# JuMinuit.jl

A native-Julia port of [CERN ROOT Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html),
the workhorse function-minimization library used throughout high-energy
physics for χ² and likelihood fits.

## Why?

[iminuit](https://github.com/scikit-hep/iminuit) (Python) and
[IMinuit.jl](https://github.com/fkguo/IMinuit.jl) (Julia, by the same
lead author) both wrap the upstream C++ library. **JuMinuit.jl is a clean-room
Julia port** of the same algorithms — no C++ dependency, no FFI overhead,
and full access to Julia tooling (ForwardDiff, threads, broadcasted FCN
evaluation). On the benchmark corpus it runs in the **0.13–0.89× C++
wall-time** range, i.e. comparable to or faster than C++ Minuit2 — see
[`benchmark/`](https://github.com/fkguo/JuMinuit.jl/tree/main/benchmark).

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

Or the iminuit / IMinuit.jl-style front end, with named parameters, limits,
and asymmetric MINOS errors:

```julia
m = Minuit(x -> sum(abs2, x .- [1.0, 2.0, 3.0, 4.0]), zeros(4);
           names = ["a", "b", "c", "d"])
migrad!(m)
minos!(m)
m.values        # ≈ [1, 2, 3, 4]
m.minos_errors  # asymmetric ±σ per parameter
```

## What's included

- **Minuit2 algorithms** — MIGRAD, HESSE, MINOS, MnContours, Simplex and
  Scan; bounds, fixed parameters, and Strategy levels 0/1/2, ported with
  line-by-line C++ fidelity and iminuit-matching defaults.
- **iminuit / IMinuit.jl-compatible front end** — `m.values`, `m.errors`,
  `migrad!`, `minos!`, `mncontour`, named-parameter access, per-parameter
  `fix!`/`set_limits!`, and Jupyter-first rich output. `Fit`/`ArrayFit` are
  exported aliases of [`Minuit`](@ref).
- **[Cost functions](cost_functions.md)** — a Julia-native family
  (`LeastSquares`, `UnbinnedNLL`, `BinnedNLL`, …) composable with `CostSum`.
- **[Error analysis](error_analysis.md) beyond HESSE/MINOS** — Monte-Carlo
  Δχ² regions, bootstrap, jackknife, and multi-modal solution detection, for
  the flat or strongly non-Gaussian likelihoods where MINOS struggles.
- **AD & threaded gradients** — a ForwardDiff extension and an opt-in
  threaded numerical gradient — plus an `Optim.jl` alternative-minimizer
  bridge (`scipy`).

## Tutorials & reference

* [Quickstart](tutorials/quickstart.md) — a hands-on tour.
* [Bounded parameters](tutorials/bounded.md) — parameter limits and fixed
  parameters.
* [MINOS errors & contours](tutorials/minos_contours.md) — asymmetric error
  bars and 2-D confidence contours.
* [Cost functions](cost_functions.md) — the Julia-native cost family.
* [Error analysis](error_analysis.md) — which uncertainty method to use, when.
* Full [API reference](api.md) and [internals](internals.md).

## Citation

If you use JuMinuit.jl in a publication, please also cite upstream
Minuit2 (which JuMinuit ports algorithmically):

> F. James, M. Roos, "Minuit: A System for Function Minimization and
> Analysis of the Parameter Errors and Correlations", Comput. Phys.
> Commun. 10 (1975) 343-367. https://doi.org/10.1016/0010-4655(75)90039-9

## License

LGPL-2.1-or-later — matches upstream Minuit2 (the same algorithms,
ported to Julia). See [`docs/UPSTREAM.md`](https://github.com/fkguo/JuMinuit.jl/blob/main/docs/UPSTREAM.md)
for provenance and attribution.
