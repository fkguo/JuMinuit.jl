# SPDX-License-Identifier: LGPL-2.1-or-later

"""
    CostFunction{F,T}

Wraps a user FCN `f::F` with an error definition `up::T` and an internal
call counter. **Closure-specialized** via parametric `F` (ROADMAP §2.3
+ Risk #4) so the call site `cf(x)` devirtualizes through the concrete
type rather than hitting Julia's `::Function` vtable.

Mirrors `MnFcn` and `MnUserFcn` from C++ Minuit2, collapsed since Julia
doesn't need the C++ inheritance hierarchy (the call counter is in the
same struct as the user function; multiple dispatch handles
gradient-vs-no-gradient via a separate `CostFunctionWithGradient` type
which lands in Phase 1).

In Phase 0 the call to `f(x)` passes the parameter vector unchanged
(identity transform: no bounds). Phase 1 inserts the internal→external
sin/sqrt transform via [`Transformation`](@ref) on the same call boundary.

# Fields

- `f::F` — the user function. Must accept an `AbstractVector{Float64}`
  (or compatible) and return a `Float64`-convertible cost value.
- `up::T` — error definition (`1.0` for χ² fits, `0.5` for negative
  log-likelihood fits; default `1.0`). Parametric `T` keeps the door
  open for `ForwardDiff.Dual{...,Float64}` users in Phase 2.1 without
  re-shuffling the type.
- `nfcn::Base.RefValue{Int}` — call counter; mutates on each
  invocation via the call-operator overload. `Ref` is the idiomatic
  Julia way to hold mutable state inside an otherwise-immutable struct.

# Performance notes

- This struct is **not** `isbits` (because of the `Ref`). One heap
  allocation per `CostFunction` instance. That's fine — we construct
  one per `migrad` call.
- The `f(x)::Float64` annotation enforces the return type contract at
  the call boundary; if the user's FCN returns a non-`Float64`, you'll
  see a clear runtime error at first call, not silent type instability
  in the MIGRAD inner loop.

# Examples

```julia
julia> cf = CostFunction(x -> sum(abs2, x), 1.0);

julia> cf([1.0, 2.0, 3.0])
14.0

julia> ncalls(cf)
1

julia> reset_ncalls!(cf); ncalls(cf)
0
```
"""
struct CostFunction{F,T}
    f::F
    up::T
    nfcn::Base.RefValue{Int}
end

CostFunction(f, up = 1.0) = CostFunction(f, up, Ref(0))

"""
    (cf::CostFunction)(x::AbstractVector) -> Float64

Evaluate the user function at `x`, incrementing the call counter.
Returns `Float64`. Numeric returns (e.g. `Int`) are coerced via the
`Float64` constructor — isbits→isbits, zero allocation. Non-numeric
returns trigger a `MethodError`.
"""
@inline function (cf::CostFunction)(x::AbstractVector)
    cf.nfcn[] += 1
    return Float64(cf.f(x))::Float64
end

"""
    ncalls(cf::CostFunction) -> Int

Number of times this `CostFunction` has been called.
"""
ncalls(cf::CostFunction) = cf.nfcn[]

"""
    reset_ncalls!(cf::CostFunction) -> CostFunction

Reset the call counter to zero. Returns `cf`.
"""
function reset_ncalls!(cf::CostFunction)
    cf.nfcn[] = 0
    return cf
end

"""
    errordef(cf::CostFunction)

Return the error definition (`up`) — `1.0` for χ², `0.5` for NLL.
"""
errordef(cf::CostFunction) = cf.up
