# Changelog

All notable changes to JuMinuit.jl. Follows [Keep a Changelog](https://keepachangelog.com/)
+ [Semantic Versioning](https://semver.org/).

## [0.1.0-alpha] — 2026-05-25

First substantial alpha release. Phase 0 PoC + Phase 1 batch 1-3 +
Phase 2.1/2.4/2.5 + Phase 3 first cut shipped. 35 commits. 888/888
tests passing. Aqua + JET clean.

### Added

#### Phase 0 — Core MIGRAD
- `MachinePrecision`, `Strategy` (levels 0/1/2 with C++-exact constants).
- `MinimumState`, `CovStatus` enum, four state types (Parameters,
  Error, Gradient, full State).
- `CostFunction{F,T}` with parametric closure specialization +
  `Ref{Int}` call counter.
- Symmetric storage convention (`:U`) + BLAS-backed kernels in
  `linalg.jl`: `sym_mul!`, `sym_rank1_update!`, `sym_invert!`,
  `sym_eigvals`, `sum_sym`, `add_sym!`.
- Numerical gradient (`InitialGradientCalculator` cold-start +
  `Numerical2PGradientCalculator` two-point central-diff refinement).
- DFP Hessian update (rank-2 base + additive rank-1 correction when
  `delgam > gvg`).
- EDM (Expected Distance to Minimum) estimator.
- `MnPosDef` positive-definiteness enforcement via eigenvalue
  perturbation.
- `NegativeG2LineSearch` for the seed when initial g2 has negative
  entries.
- Parabolic 1D line search (`MnLineSearch` minus the deferred
  cubic/Brent variants).
- `MnSeedGenerator` for the initial MIGRAD state.
- Full MIGRAD loop with `FunctionMinimum` result type.

#### Phase 1 — Bounds + MINOS + Contours + HESSE
- `transform.jl`: sin / SqrtUp / SqrtLow / identity bound transforms
  matching C++ exactly (including the sign-aware SqrtUp derivative).
- `parameters.jl`: `MinuitParameter` + `Parameters` (collapsed
  `MnUserParameters` + `MnUserTransformation`).
- `hesse.jl`: full numerical Hessian with diagonal multiplier loop,
  off-diagonal pass, `MnPosDef` + invert, status flag handling.
- `covariance_squeeze.jl`: drop a row+col from a symmetric matrix via
  invert → squeeze → invert back, with diagonal fallback on failure.
- `function_cross.jl`: parabolic root-find with inner re-minimization;
  used by MINOS.
- `minos.jl`: asymmetric ±σ errors with `MinosError` result type.
- `contours.jl`: 2D 1σ contour via ellipse approximation from MINOS +
  off-diagonal covariance.
- `migrad_bounded.jl`: bound-aware MIGRAD via `Parameters` wrapper;
  internal MIGRAD operates in unbounded coords, user FCN sees
  external coords; full external covariance back-conversion via
  Jacobian chain rule.

#### Phase 2 — Polish
- `ad_gradient.jl` (2.1): `CostFunctionWithGradient{F,G,T}` for
  user-supplied or AD-produced gradients; ForwardDiff integration.
- `serialize.jl` (2.5): `to_dict` / `minimum_summary_from_dict` for
  JSON / JLD2 roundtrip of all result types.
- `precompile_workload.jl` (2.4): PrecompileTools workload reducing
  TTFX by ~50% on typical MIGRAD paths.

#### Phase 3 — User API
- `minuit.jl`: iminuit-style `Minuit` mutable struct with
  `migrad!`, `minos!`, `contour` methods and `m.values`, `m.errors`,
  `m.fval`, `m.edm`, `m.nfcn`, `m.valid`, `m.covariance` property
  access (via `Base.getproperty`).

#### Tooling
- `tools/cpp_trace_harness.cxx`: C++ Minuit2 reference-data generator
  producing JSON oracles for unbounded + bounded + fixed-parameter
  benchmark cases.
- `tools/regen_reference.sh`: build + run wrapper.
- `benchmark/cpp/cpp_bench.cxx`: wall-time benchmark of C++ Minuit2
  for §3.4 Criterion 2 cross-implementation comparison.
- `benchmark/compare_cpp.jl`: pulls Julia + C++ medians, prints
  ratio table, computes verdict.
- `benchmark/bench_migrad_suite.jl` + `benchmark/perf-config.toml`:
  julia-perf Level-2 evidence-gate suite.
- `scripts/run_gate.sh`: gate driver.

### Verified

- **Phase 0 §3.4 Criterion 1**: Quad-4D matches C++ Minuit2 reference
  JSON to fval ≤ 1e-15, params to 1e-10. Rosenbrock cases within
  Strategy(0) cross-impl variance.
- **Phase 0 §3.4 Criterion 2**: Julia ≤ 0.887× C++ wall time on
  every benchmark in the §3.3 corpus (max ratio 0.887×, mean 0.47×).
- **Phase 0 §3.4 Criterion 4**: Aqua + JET clean on the public API
  (`migrad(::Function, ::Vector{Float64}, ::Vector{Float64})`).
- **Phase 1 bounded oracle parity**: 4 bounded reference cases (Sin /
  upper-only / lower-only / fixed-parameter) match C++ Minuit2 output
  on fval, free-parameter values, and NFcn within documented Strategy(0)
  tolerance. External covariance verified symmetric.

### Audit trail

- Four rounds of independent parallel multi-agent review (codex
  gpt-5.5 xhigh + native Opus subagent), all archived under
  `scratch/{codex,opus}_review_phase*.md`:
  1. v1 → v2 ROADMAP reconciliation (caught a real `sum_sym`
     signed-vs-absolute blocking bug in linalg).
  2. Phase 0 MIGRAD integration (caught `reached_call_limit` AND-gate
     bug, plus 9 surgical issues).
  3. Phase 0 hot-path kernel retroactive review.
  4. Phase 1 mid-phase (caught the `initial_int_errors` Taylor vs
     two-sided perturbation divergence at bounds, plus 5 minor).
  5. Phase 1 batch 2+3 (caught D3 covariance asymmetric-read bug, C-2
     contour sign-blind selector, A7/B4 bounded coord-frame leak,
     A5 strategy no-op, B2 MinosError field semantics).

  All blocking and high-priority findings applied as commits with
  source-cited diffs.

### Deferred

#### Phase 1.x
- `function_cross` C++ 3-point parabolic algorithm parity (Julia
  uses a simplified 2-point fit + replace-worst).
- Multi-parameter `function_cross` for the C++-exact (non-ellipse)
  contour algorithm.
- Strategy(1+) `HessianGradientCalculator` refinement inside
  `MnHesse`.
- Bounded MINOS through the internal-coord-wrapped CostFunction
  (currently the unbounded MINOS path is wired; bounded MINOS via
  `Minuit.minos!(m, par)` uses the wrapped CF but doesn't refine
  bounds at the parabolic-cross step).
- `Int2extError` two-sided bounded errors (currently uses
  Jacobian-diagonal sqrt; matters near bounds).
- Variable-sized external covariance to match C++
  `MnUserParameterState` shape (currently full n_total × n_total
  with zero rows for fixed parameters).

#### Phase 2
- 2.2 Threads-parallel numerical gradient.
- 2.3 Plot recipes (RecipesBase).

#### Phase 3 polish
- Full iminuit pretty-print parity.
- `m.errors[name] = ...` setter API.
- Documentation site (Documenter.jl).

[0.1.0-alpha]: https://github.com/fkguo/JuMinuit.jl/compare/main...HEAD
