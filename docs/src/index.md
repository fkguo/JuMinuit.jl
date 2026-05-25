# JuMinuit.jl

A native-Julia port of [CERN ROOT Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html),
the workhorse function-minimization library used throughout high-energy
physics for χ² and likelihood fits.

## Why?

[iminuit](https://github.com/scikit-hep/iminuit) (Python) and
[IMinuit.jl](https://github.com/fkguo/IMinuit.jl) (Julia, by the same
author) both wrap the upstream C++ library. **JuMinuit.jl is a clean-room
Julia port** of the same algorithms — no C++ dependency, no FFI overhead,
and full access to Julia tooling (ForwardDiff, threads, broadcasted FCN
evaluation).

Performance target: **comparable to or better than C++** on the Phase 0
benchmark corpus (Rosenbrock, Quad-NF, Gauss-LL). On Apple Silicon with
4 threads the current code is in the 0.13–0.89× C++ wall-time range —
see [`benchmark/`](https://github.com/fkguo/JuMinuit.jl/tree/main/benchmark).

## Quick example

```julia
using JuMinuit

# χ² with a 4-parameter quadratic
cf = CostFunction(x -> sum(abs2, x .- [1.0, 2.0, 3.0, 4.0]))

# Initial parameter values + step sizes
m = migrad(cf, [0.0, 0.0, 0.0, 0.0], [0.1, 0.1, 0.1, 0.1])

show(stdout, MIME"text/plain"(), m)
```

```
┌───────────────────────────────────────────────────────────────────────┐
│                                Migrad                                 │
├───────────────────────────────────┬───────────────────────────────────┤
│ FCN = 4.196e-16                   │             Nfcn = 18             │
│ EDM = 4.196e-16 (Goal: 0.002)     │                                   │
├───────────────────────────────────┼───────────────────────────────────┤
│           Valid Minimum           │  Below EDM threshold (goal x 10)  │
├───────────────────────────────────┼───────────────────────────────────┤
│      No parameters at limit       │         Below call limit          │
├───────────────────────────────────┼───────────────────────────────────┤
│             Hesse OK              │        Covariance accurate        │
└───────────────────────────────────┴───────────────────────────────────┘
```

## Status

| Phase | What it ships | Status |
|------:|:--------------|:-------|
| 0 | MIGRAD (Strategy 0), gradient calculator, DFP update, EDM, linesearch, posdef | ✅ done |
| 1 | Bounds (sin/sqrt transforms), parameters API, HESSE, MINOS, contours, covariance squeeze, inner-Hesse for Strategy ≥ 1 | ✅ done |
| 1.x | C++-exact MnFunctionCross 3-point parabolic, Int2extError two-sided, free-covariance accessor | ✅ done |
| 2.1 | AD-backed analytical gradients (ForwardDiff) | ✅ done |
| 2.2 | Threaded numerical gradient | ✅ done |
| 2.3 | Plots/RecipesBase recipes | ✅ done |
| 2.4 | PrecompileTools workload | ✅ done |
| 2.5 | JSON serialization | ✅ done |
| 3 | iminuit-style pretty-print, Documenter docs site, polish | ⏳ in progress |

## Next steps

* See [Quickstart](tutorials/quickstart.md) for a hands-on tour.
* For parameter limits and fixed parameters, see [Bounded parameters](tutorials/bounded.md).
* For asymmetric error bars and 2D confidence contours,
  see [MINOS errors & contours](tutorials/minos_contours.md).
* Full [API reference](api.md) and [internals](internals.md).

## Citation

If you use JuMinuit.jl in a publication, please also cite upstream
Minuit2 (which JuMinuit ports algorithmically):

> F. James, M. Roos, "Minuit: A System for Function Minimization and
> Analysis of the Parameter Errors and Correlations", Comput. Phys.
> Commun. 10 (1975) 343-367. https://doi.org/10.1016/0010-4655(75)90039-9

## License

LGPL-2.1-or-later — matches upstream Minuit2 (the same algorithms,
ported to Julia).
