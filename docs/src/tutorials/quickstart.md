# Quickstart

This tutorial walks through a minimal MIGRAD fit and the key result
accessors. We assume Julia ≥ 1.10 and `JuMinuit` already installed
(`Pkg.add("JuMinuit")` when registered, or `Pkg.add(url=...)` from
the repository).

## A χ² fit in five lines

```julia
using JuMinuit

cf = CostFunction(x -> sum(abs2, x .- [1.0, 2.0, 3.0]))   # χ² = (x-1)² + (y-2)² + (z-3)²
m  = migrad(cf, [0.0, 0.0, 0.0], [0.1, 0.1, 0.1])         # initial point + step sizes
@assert is_valid(m)
println("fval = ", fval(m))                                # 0 at the minimum
println("x*   = ", Base.values(m))                         # [1, 2, 3]
```

## Anatomy of `CostFunction`

`CostFunction(f, up=1.0)` wraps a Julia function `f(x::AbstractVector{<:Real}) -> Real`
and an "error definition" `up`:

| `up`  | When to use                                          |
|------:|:-----------------------------------------------------|
| `1.0` | χ² fits (default). 1σ contour at `f = fmin + 1`.     |
| `0.5` | Negative-log-likelihood fits. 1σ at `f = fmin + 0.5`.|

The `up` value matters for MINOS and contour geometry, not for the
location of the minimum.

## Accessing results

After a successful `migrad`, the returned `FunctionMinimum` exposes:

```julia
fval(m)         # value of the FCN at the minimum
Base.values(m)  # Vector{Float64} of parameter values at the minimum
errors(m)       # Vector of 1σ Hesse errors (sqrt(2·up·V[i,i]))
covariance(m)   # Symmetric{Float64} covariance matrix (2·up·V)
gradient(m)     # FunctionGradient at the minimum
nfcn(m)         # total FCN calls used by MIGRAD
is_valid(m)     # true if MIGRAD converged within tolerances
```

## Tolerances and budgets

```julia
m = migrad(cf, x0, errs;
    tol = 0.1,           # EDM convergence multiplier (× up × 0.002)
    maxfcn = 1000,       # call budget
    strategy = Strategy(0),     # 0 = fast, 1 = default, 2 = thorough
    prec = MachinePrecision(),  # numerical precision overrides
)
```

`Strategy(0)` is the speed-tuned mode (no inner-Hesse refinement).
`Strategy(1)` (the iminuit default) runs an inner HESSE when the
DFP-estimated covariance is loose (`Dcovar > 0.05`). `Strategy(2)`
runs HESSE unconditionally after MIGRAD converges and re-iterates the
inner DFP loop if HESSE moves EDM above tolerance.

## Convergence diagnostics

If `is_valid(m)` returns `false`, inspect:

```julia
m.reached_call_limit  # true if maxfcn exhausted
m.above_max_edm       # true if EDM > 10 × tol·up·0.002
m.hesse_failed        # true if Strategy ≥ 1 inner HESSE failed
m.made_pos_def        # true if posdef perturbation was applied mid-fit
```

For invalid results, try (in order):

1. Loosen `tol`.
2. Increase `maxfcn`.
3. Raise the strategy level.
4. Re-seed with a better starting point.

## Next

Continue to [Bounded parameters](bounded.md) to learn how to add
`lower`/`upper` limits and `fixed` parameters.
