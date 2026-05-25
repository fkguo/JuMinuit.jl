# Bounded parameters

Phase 1 ships the iminuit-style **named-parameter API**: `MinuitParameter`
+ `Parameters`, with per-parameter `lower`/`upper`/`fixed` flags.
Bounds are implemented via internal `sin`/`sqrt` coordinate transforms
mirroring C++ Minuit2 (no L-BFGS-B-style projected gradient).

## Adding bounds

```julia
using JuMinuit

cf = CostFunction(x -> (x[1] - 0.5)^2 + (x[2] - 0.5)^2)

params = Parameters([
    MinuitParameter("x", 0.3, 0.1; lower = 0.0, upper = 1.0),  # double-bounded
    MinuitParameter("y", 0.3, 0.1),                              # unbounded
])

m = migrad(cf, params)
@assert is_valid(m)
@assert 0.0 ≤ m.ext_values[1] ≤ 1.0   # bound respected exactly
```

The first positional argument to `MinuitParameter` is the **name**,
followed by the initial value and step size. Keywords:

| Keyword     | Default | Meaning                                                  |
|:------------|:--------|:---------------------------------------------------------|
| `lower`     | `NaN`   | Lower bound; absent if `NaN`.                            |
| `upper`     | `NaN`   | Upper bound; absent if `NaN`.                            |
| `fixed`     | `false` | If `true`, parameter stays at its initial value.         |

When both `lower` and `upper` are set, the parameter uses the **Sin**
transform: `ext = lower + 0.5·(upper-lower)·(sin(int)+1)`. With only
one bound, the SqrtUp or SqrtLow transform is used. The internal
coordinate the optimizer sees is **unbounded** — projection onto the
allowed range happens transparently in `int2ext`/`ext2int`.

## Fixed parameters

```julia
params = Parameters([
    MinuitParameter("x", 0.0, 0.1),
    MinuitParameter("y", 5.0, 0.1; fixed = true),    # bit-exact at 5.0
    MinuitParameter("z", 0.0, 0.1),
])

m = migrad(cf, params)
@assert m.ext_values[2] == 5.0          # exact, no roundoff
@assert m.ext_errors[2] == 0.0           # fixed → zero error
```

## Result accessors for the bounded path

```julia
m = migrad(cf, params)         # returns BoundedFunctionMinimum

m.ext_values                   # parameter values in external (user) coords
m.ext_errors                   # symmetric-average two-sided errors via Int2extError
ext_covariance(m)              # full n_total × n_total external covariance
                               #   (fixed rows + cols are zero)
free_covariance(m)             # n_free × n_free sub-block (C++ MnUserParameterState shape)
```

The diagonal of `m.ext_errors` uses the C++
[`MnUserTransformation::Int2extError`](https://github.com/root-project/root/blob/master/math/minuit2/src/MnUserTransformation.cxx)
two-sided formula:

```
ui  = int2ext(val)
du1 = int2ext(val + err) − ui
du2 = int2ext(val − err) − ui
return 0.5 · (|du1| + |du2|)
```

This captures the nonlinear remapping near bounds where the
Jacobian-only `sqrt(V_ext[i,i])` would under-report (the Jacobian
shrinks toward zero at the boundary).

## Bound saturation

If MIGRAD ends with a parameter pinned at its `lower` or `upper` limit,
the pretty-print will reflect "Some parameters at limit" (Phase 3
extension; current minimal output flags this in the box header).
