# JuMinuit.jl Roadmap

A native-Julia port of the C++ Minuit2 function minimization library
(GooFit/Minuit2; ROOT::Math::Minuit2), targeting drop-in replacement of
the iminuit/IMinuit.jl stack with C++-comparable performance.

> **Compass**: this document is read side-by-side with the C++ reference
> at `reference/Minuit2_cpp/`. Filenames cited without a path live there.
> The module-mapping table (§7) is the authoritative porting compass.

---

## 1. Project goal

Reimplement, in idiomatic Julia, the Variable-Metric MIGRAD algorithm
and its companions (SIMPLEX, MINOS, HESSE, CONTOURS) that have defined
HEP fitting for forty years, while matching the C++ reference implementation
to ~1e-10 on shared benchmarks. The port eliminates the C++ runtime dependency
that IMinuit.jl currently carries through PyCall/iminuit, gives Julia-native
users automatic differentiation, threaded likelihoods, and Plots/Makie recipes,
and keeps the user-visible API close to iminuit's so existing fits migrate.

The Variable-Metric MIGRAD loop is a thirty-year-old optimized inner cycle.
We will **port it, not redesign it**. Numerical reproducibility against the C++
reference is the primary acceptance criterion; performance follows.

---

## 2. Performance philosophy

Minuit2's hot paths are very stereotyped. Knowing where the cycles go in the C++
code tells us exactly where Julia must be careful.

### 2.1 Where the C++ spends its time

For a fit with `n` free parameters per MIGRAD iteration:

| Hot path                       | C++ location                              | C++ cost / iter |
|--------------------------------|-------------------------------------------|-----------------|
| User FCN evaluations           | `MnUserFcn::operator()` → user lambda     | ~2n + line-search calls (5–15) |
| Numerical gradient (central diff, multi-cycle) | `Numerical2PGradientCalculator.cxx:63–230` | `2·n·Ncycle` FCN calls — **dominant** for cheap FCNs |
| Step = -V·g (DSPMV)            | `LAVector.h:122–131`, `mndspmv.cxx`       | O(n²) FLOPs, packed sym |
| Inner products `g·g`, `step·g` | `LaInnerProduct.cxx`, `mnddot.cxx`        | O(n) |
| DFP Hessian update             | `DavidonErrorUpdator.cxx:24–73`           | O(n²) — outer products, scalar adds |
| EDM = ½ gᵀ V g                 | `VariableMetricEDMEstimator.cxx`          | O(n²) (similarity) |
| Line search                    | `MnLineSearch.cxx`                        | 4–12 FCN calls, scalar math |
| Final HESSE (numerical 2nd derivs) | `MnHesse.cxx`                          | `O(Ncycle·n) + n(n-1)/2` FCN calls — one-shot |

Three regimes determine the optimization target:

- **Cheap FCN** (e.g. Rosenbrock, low-dim quadratic): gradient and linear-algebra
  loops dominate — Julia must match C++ on dense Float64 BLAS and zero-allocate
  the inner loop.
- **Moderate FCN** (1k–10k events): FCN cost ≈ gradient overhead — type stability
  of the user closure is critical.
- **Expensive FCN** (large unbinned likelihoods): user code dominates — Julia
  shines (closures specialize), and Threads-parallel reductions over data may
  win outright vs. single-threaded C++.

### 2.2 Concrete Julia idioms (and what they map to in C++)

**Type stability & inference**
- `FCN` is held as a parametric type: `struct CostFunction{F,T} ; f::F ; up::T ; end`
  — never as `Function` (the C++ `FCNBase &` is a vtable call; Julia closure call
  through a concrete type devirtualizes).
- All state structs (`MinimumState`, `FunctionGradient`, `MinimumError`,
  `MinimumParameters`) are concrete, parameterized only on element type when
  needed. No `Any` fields.

**In-place linear algebra (zero allocation in the MIGRAD loop)**
- The MIGRAD loop in `VariableMetricBuilder.cxx:237–341` allocates two
  per-iteration temporaries (`step`, `vUpd`) and one DFP update; everything else
  is in-place via BLAS calls. Julia must mirror:
  - Preallocate `step::Vector{Float64}`, `g::Vector{Float64}`, `g_prev::Vector{Float64}`,
    `dx::Vector{Float64}`, `dg::Vector{Float64}`, `vg::Vector{Float64}` once in
    a `MigradWorkspace` struct attached to the builder.
  - Replace `vUpd = Outer_product(dx)/delgam - Outer_product(vg)/gvg` with
    `BLAS.spr!` / hand-rolled symmetric outer-product into a preallocated
    `Vector{Float64}` packed-storage buffer; OR, use a plain `Matrix{Float64}`
    + `Symmetric` view and call `BLAS.syr!` (LAPACK preference) — choose based
    on benchmarks (see §3 exit criterion).
  - `step = -V * g` → `BLAS.spmv!('U', -1.0, V_packed, g, 0.0, step)` or
    `mul!(step, Symmetric(V), g, -1.0, 0.0)`.
- `@views` for any slice; `@inbounds` after correctness tests pass; avoid
  building intermediate `Vector`s in arithmetic.

**Small-dimension specialization**
- Profile-driven: if a meaningful fraction of usage is n ≤ ~16, we may add an
  `SVector`/`SMatrix` (StaticArrays.jl) inner path. **Defer** — measure first.
  C++ does not specialize for small n; only StackAllocator helps it there.

**Memory management**
- StackAllocator (`StackAllocator.h`) is irrelevant in Julia — its purpose is to
  bypass `malloc`. Julia's GC pressure is solved by the preallocation strategy
  above, not by a custom allocator.
- All scratch buffers live in workspace structs passed through the call chain;
  no globals.

**LAPACK directly when beneficial**
- `mnvert.cxx` is a Gauss–Jordan symmetric inversion. For `n ≥ ~8`, switching to
  `LAPACK.sptrf!` / `LAPACK.sptri!` (Bunch–Kaufman on packed symmetric) is both
  faster and numerically better. Reasonable threshold to be measured.
- `LASymMatrix` uses lower-triangular packed storage; LAPACK's `*sp*` family
  is the direct equivalent.

**Eigenvalue routines (MnPosDef)**
- `eigenvalues(p)` in `MnPosDef.cxx:80` → `LAPACK.spev!('N','U', packed)` or
  `eigvals(Symmetric(M))`.

**Parallel FCN evaluation**
- The C++ uses OpenMP in `Numerical2PGradientCalculator.cxx:112–127` to
  parallelize the gradient loop. Julia: `Threads.@threads` over the parameter
  index, **but only as a Phase 2 feature** — Phase 0 stays single-threaded for
  apples-to-apples comparison.

**Type-stable error definitions**
- `up::Float64` lives in the FCN struct; we don't carry it through the call
  chain dynamically.

### 2.3 Pitfalls specific to this port

- **Closure specialization**: a `FCN(::AbstractVector) -> Real` user closure must
  be held via parametric `F`. A non-parametric `Function` field will silently
  bring vtable-grade overhead.
- **`@views` vs new alloc on slicing**: `g[k] - g_prev[k]` works element-wise but
  `g - g_prev` allocates. Use `@. dg = g - g_prev`.
- **BLAS thread interaction**: when the user wires multi-threaded likelihoods
  (Phase 2), `BLAS.set_num_threads(1)` is the safe default to avoid nested
  parallelism — same trick as iminuit + numpy.
- **C++ ABObj expression templates** (`ABObj.h`, `ABSum.h`, `ABProd.h`): exist
  *only* to fuse `step = a*M*v + b*w` into one BLAS call without temporaries.
  Julia's broadcasting + `mul!` / `axpy!` solve this cleanly without the template
  machinery. **Do not port ABObj**; replace with explicit `mul!`/`axpy!` calls.
- **`shared_ptr<BasicX>` indirection**: `MinimumState`, `FunctionGradient`, etc.
  use shared-pointer handles in C++ (e.g. `FunctionMinimum.h:99`). In Julia, plain
  immutable `struct`s with `mutable` only where genuine state mutation matters
  (the workspace). Don't blindly mirror the shared-handle layer.

---

## 3. Phase 0 — Proof of concept

**Mandate**: prove the Julia port can reach C++-Minuit2 performance on simple
MIGRAD-only fits, within numerical equivalence, before any large code surface is
written.

**Scope (in)**: unconstrained MIGRAD (gradient-based) with numerical gradient,
free parameters only, no bounds, no fixed, no MINOS, no SIMPLEX. Enough to fit
Rosenbrock and a Gaussian negative-log-likelihood.

**Scope (out)**: parameter bounds (sin/sqrt transformations), fixed parameters,
MINOS errors, contours, SIMPLEX, Fumili, the IMinuit.jl API surface, anything
plot-related.

**Exit gate**: ≤ 1.5× C++ wall time on the Rosenbrock-10 and Gauss-LL-100
benchmarks, results identical to 1e-10 on parameter values and function minimum.
If the gate fails, the architectural assumptions in §2 are wrong and must be
revisited before Phase 1 starts.

### 3.1 Files to create under `src/`

src/
  JuMinuit.jl                # top-level module; reexports the public surface
  precision.jl               # MachinePrecision: mirrors MnMachinePrecision
  strategy.jl                # Strategy: mirrors MnStrategy (levels 0/1/2)
  fcn.jl                     # CostFunction{F} wrapper; ncalls counter
  parameters.jl              # MinuitParameter (single param); free-only first
  state.jl                   # MinimumParameters, MinimumError, FunctionGradient,
                             # MinimumState; concrete, immutable
  workspace.jl               # MigradWorkspace: all scratch buffers
  linalg.jl                  # symmetric packed↔dense helpers, mul!, axpy! wrappers,
                             # invert!, eigvals — thin layer for benchmarking
  gradient.jl                # initial_gradient!, numerical_gradient! (Numerical2P)
  edm.jl                     # estimate_edm (mirrors VariableMetricEDMEstimator)
  posdef.jl                  # make_posdef! (mirrors MnPosDef)
  davidon.jl                 # davidon_update! (DFP, with rank-1 branch)
  linesearch.jl              # parabolic line search (mirrors MnLineSearch op())
  migrad.jl                  # the iteration loop (mirrors VariableMetricBuilder)
  seed.jl                    # seed_state (mirrors MnSeedGenerator)
  result.jl                  # FunctionMinimum
  api.jl                     # migrad(fcn, x0, errors; strategy, tol, maxfcn)
```

**Why this many files**: each file maps 1-to-1 to a C++ translation unit we'll
diff against during the port. Smaller modules make line-by-line review tractable
and isolate any numerical regression to a unit test.

**What NOT to create yet**:
- `transform.jl` (sin/sqrt bound transforms) — Phase 1
- `migrad_api.jl` mimicking `MnMigrad` — Phase 1
- `minos.jl`, `hesse.jl`, `contours.jl` — Phase 1 (HESSE) / Phase 1 (MINOS)
- `simplex.jl` — Phase 1
- `precompile.jl`, plotting, AD glue — Phase 2

### 3.2 Tests under `test/`

TDD with the C++ reference as oracle. Test files mirror sources:

```
test/
  runtests.jl                # top-level driver
  test_linalg.jl             # spmv, syr, invert; cross-checked vs LAPACK
  test_precision.jl          # MachinePrecision constants match C++
  test_strategy.jl           # strategy level 0/1/2 values match MnStrategy.cxx:33-70
  test_initial_gradient.jl   # InitialGradientCalculator outputs on Quad4
  test_numerical_gradient.jl # Numerical2P outputs on Rosenbrock & Quad4;
                             # exact reproduction at chosen seed
  test_davidon.jl            # synthetic 4D Hessian update reproduces C++ to 1e-12
  test_linesearch.jl         # canned step/gradient configurations
  test_edm.jl                # gᵀ V g / 2 numerical equivalence
  test_migrad_quad4.jl       # full MIGRAD on Quad4F (MnTutorial/Quad4F.h);
                             # asserts min == 0, params == 0 to 1e-10, same NFcn
  test_migrad_rosenbrock.jl  # 2D and 10D Rosenbrock; results vs reference dump
  test_migrad_gauss_ll.jl    # 100-point Gaussian NLL fit (mirror of
                             # MnSim/GaussFcn.cxx); compare to fixed-seed reference
  reference_data/            # JSON dumps of C++ runs (one-time generation)
    quad4f_min.json
    rosenbrock2d_min.json
    rosenbrock10d_min.json
    gaussll_min.json
```

**Reference data generation** (one-time, by hand or in `tools/`):
Build the C++ examples under `reference/Minuit2_cpp/examples/simple/` (the
`Quad1F` example shows the build pattern). Add a small CMakeLists wrapper
that prints `(min.Fval(), min.Edm(), min.NFcn(), each param value/error)` to
JSON. Check the JSON into `test/reference_data/` so the test suite stays
self-contained and doesn't require a working C++ toolchain.

### 3.3 Benchmark suite under `benchmark/`

```
benchmark/
  bench_migrad.jl            # BenchmarkTools.jl; Rosenbrock, Quad4, Gauss-LL
  bench_gradient.jl          # isolated numerical gradient on Rosenbrock-10
  bench_davidon.jl           # isolated DFP update at varying n
  cpp/                       # tiny CMake project that builds the same fits
    CMakeLists.txt
    bench_rosenbrock.cxx
    bench_gauss_ll.cxx
  compare.jl                 # runs Julia BenchmarkTools median, parses C++ stderr
                             # wall-time output, prints a ratio table
  results/                   # tracked: history of (commit, machine, ratio) JSON
```

**Benchmark scenarios** (the exit-gate corpus):

| Name                  | n free | FCN cost   | Why it matters |
|-----------------------|--------|------------|----------------|
| Rosenbrock-2          | 2      | trivial    | low-n, gradient-loop dominated |
| Rosenbrock-10         | 10     | trivial    | the canonical Minuit stress test; non-quadratic curvature |
| Quad4F (analytic mode)| 4      | trivial    | sanity check, exact convergence |
| Gauss-LL-2 × 100      | 2      | 100 events | typical "small fit", FCN ≈ overhead |
| Gauss-LL-10 × 1000    | 10     | 1k events  | typical HEP fit |
| Gauss-LL-40 × 1000    | 40     | 1k events  | parallel-test analog (`MnSim/ParallelTest.cxx`) |

### 3.4 Acceptance criteria (Phase 0 exit gate)

A merge to `main` enabling Phase 1 requires *all* of:

1. **Correctness**: every benchmark in §3.3 reproduces the C++ reference to:
   - `|Δ fval| / max(1, |fval|) ≤ 1e-10`
   - `|Δ params_i| / max(1, |params_i|) ≤ 1e-10`
   - `|Δ edm| / edm ≤ 1e-6` (edm is loosely defined — see VariableMetricBuilder.cxx:66 where C++ multiplies by 0.002)
   - NFcn within ±2 of C++ (rounding-difference can swing one extra line-search step)
2. **Performance**: median Julia wall time / median C++ wall time ≤ 1.5 on every
   benchmark in §3.3, on the developer's reference machine, with single-threaded
   BLAS on both sides.
3. **Zero allocations in the inner loop**: `@allocated` for one MIGRAD iteration
   on Rosenbrock-10 returns 0 (after compile/warmup).
4. **Clean run**: no warnings, no failed asserts, `Aqua.jl` clean,
   `JET.@report_opt` clean on the public API.

### 3.5 Implementation order

Within Phase 0, the order minimizes integration risk:

1. **Day 1–3**: scaffolding (`Project.toml`, `JuMinuit.jl`, `precision.jl`,
   `strategy.jl`, `state.jl`) + `test_precision.jl` / `test_strategy.jl` green.
2. **Day 4–7**: `linalg.jl` + `test_linalg.jl` — pick packed vs. dense
   representation here, freeze the decision.
3. **Day 8–12**: `gradient.jl` (initial + numerical) + tests on Quad4 and
   Rosenbrock matching C++ to ~1e-12.
4. **Day 13–15**: `davidon.jl`, `edm.jl`, `posdef.jl` + tests.
5. **Day 16–20**: `linesearch.jl` + tests, then `seed.jl`.
6. **Day 21–25**: `migrad.jl` — the loop. End-to-end test on Quad4F first
   (closed form), then Rosenbrock-2.
7. **Day 26–28**: full benchmark sweep, profile, optimize. Run exit gate.

Any slip past day 28 triggers a design review of §2.2 idioms.

---

## 4. Phase 1 — Core port (goal + exit only)

**Goal**: feature-complete equivalence with `MnMigrad + MnHesse + MnMinos +
MnSimplex + MnContours` for the cases iminuit covers. Bounds, fixed parameters,
named parameters, MINOS asymmetric errors, contours, HESSE-after-MIGRAD,
strategy levels 0/1/2 all working.

**Exit criteria**:
- All MnTutorial examples (`Quad1F`, `Quad4F`, `Quad8F`, `Quad12F`) reproduce.
- The `MnSim/GaussFcn`-class likelihood fits reproduce.
- Parameter bounds via `sin`/`sqrt` transforms produce identical internal
  parameter values to C++ at every iteration of a chosen reference fit
  (recorded with `MnTraceObject` equivalent).
- MINOS errors agree to 1e-8 with the C++ reference on Quad4F and a
  bounded Gauss fit.
- Performance ratio ≤ 1.3× C++ on all Phase 0 benchmarks plus a bounded
  variant of Gauss-LL-10.
- No public name change between Phase 0 and Phase 1; Phase 0 users keep working.

Phase 1 task breakdown is **deliberately deferred** until Phase 0 exits.
What we learn from achieving 1.5× will reshape the Phase 1 plan.

---

## 5. Phase 2 — Julia-native extras

Each extra ships as an independent feature, gated behind its own milestone.
Order is by user demand (the IMinuit.jl user base) rather than dependency.

- **2.1 AD-backed analytical gradients**
  - `FCNGradAdapter`-equivalent that takes a user FCN + an AD backend choice
    (`ForwardDiff` for small n, `Enzyme` for large n / hot inner loops, optional).
  - Auto-fallback to numerical gradient if AD fails.
  - Verification: AD-gradient MIGRAD on Rosenbrock-10 should reduce NFcn by
    ~`2·n·NCycle` per iteration and beat numerical gradient by ≥ 2× wall.
- **2.2 Threads-parallel numerical gradient**
  - `@threads :static` over the parameter index in `numerical_gradient!`
    (analog of the OpenMP block at `Numerical2PGradientCalculator.cxx:112–127`).
  - Workspace per thread; final reduction into shared `g`, `g2`, `gstep`.
  - Default off; opt-in via `migrad(fcn, ...; threaded_grad=true)`. Document
    BLAS thread interaction.
- **2.3 Plot recipes**
  - RecipesBase recipes for: profile likelihood (1D), contour, MINOS error
    bars; mirror `MnContours` output shape so plotting just consumes the result
    struct.
- **2.4 Precompilation & startup**
  - `PrecompileTools.@compile_workload` covering `migrad(::CostFunction{F},
    ::Vector{Float64}, ::Vector{Float64})` for a few common `F` patterns;
    measure TTFX vs. without.
- **2.5 Result serialization**
  - `FunctionMinimum` → `Dict` and `Dict` → `FunctionMinimum` for easy JSON/JLD2
    persistence; useful for the JuMinuit-vs-Minuit regression CI.

---

## 6. Phase 3 — API parity with IMinuit.jl / iminuit

The user-facing entry point should let a user copy-paste from iminuit (Python)
or IMinuit.jl with at most renaming. Parity targets:

- `Minuit(fcn, x0; name=names, error=errors, limit=limits, fix=fixed, ...)`
  constructor identical to iminuit's.
- `m.migrad()`, `m.hesse()`, `m.minos()`, `m.contour(i, j)` methods.
- Property access on the result: `m.values`, `m.errors`, `m.fmin`, `m.valid`,
  `m.covariance`, `m.params` (using `Base.getproperty` overload).
- Pretty printing via `show(::IO, ::MIME"text/plain", ::Minuit)` matching
  iminuit's table-style output line by line where possible.
- IMinuit.jl compatibility shim: optional `IMinuitCompat` submodule re-exports
  IMinuit.jl's `Migrad`/`Minos`/etc. signatures so existing user code requires
  changing only the `using` line.
- Documentation: every iminuit tutorial reproduced as a Julia example in
  `docs/src/`.

**Exit criteria**: ten randomly chosen IMinuit.jl scripts in real-world HEP fits
run, with at most a `using JuMinuit` substitution, to numerical equivalence
with the original (≤ 1e-8 on parameter values).

---

## 7. Module mapping table

The porting compass. Read row-by-row against the C++ source.

| C++ class / file | Julia module / type | Phase | Notes |
|---|---|---|---|
| `Minuit2Minimizer.h/.cxx` | `JuMinuit.Minuit` (struct) + `migrad!/hesse!/minos!` methods | 3 | The ROOT-API-style facade; only Phase 3 needs it. |
| `MnApplication.h/.cxx` | `JuMinuit.Application` (internal) | 1 | Bundles FCN + state + strategy; the Phase 1 user-facing entry. |
| `MnMigrad.h` | `migrad(fcn, x0, errors; ...)` free function (P0) and `Migrad(fcn, params; ...)` struct (P1) | 0 / 1 | Many overloads in C++; one keyword-driven Julia function. |
| `MnSimplex.h` + `SimplexMinimizer.h` + `SimplexBuilder.cxx` + `SimplexParameters.cxx` + `SimplexSeedGenerator.cxx` | `JuMinuit.simplex(fcn, x0; ...)` + `SimplexBuilder` | 1 | Nelder–Mead, no derivatives. The C++ code (SimplexBuilder.cxx:24–) is ~200 lines, idiomatic Julia translation is straightforward. |
| `MnMinos.h/.cxx` + `MnFunctionCross.cxx` + `MnCross.h` + `MinosError.h` | `JuMinuit.minos(fmin, fcn, ipar; ...)` + `MinosError` struct | 1 | Asymmetric errors via `MnFunctionCross`. |
| `MnContours.h/.cxx` + `ContoursError.h` | `JuMinuit.contour(fmin, fcn, i, j; npoints=20)` | 1 | Sits on MINOS + line search. |
| `MnHesse.h/.cxx` | `JuMinuit.hesse!(state; ...)` | 1 | Full numerical Hessian; line 100–315 of MnHesse.cxx is the dense path. |
| `MnScan.h/.cxx` + `ScanBuilder.cxx` + `MnParameterScan.cxx` | `JuMinuit.scan(fcn, ipar; ...)` | 2 | 1D function scan; mostly cosmetic. |
| `ModularFunctionMinimizer.h/.cxx` | (folded into `migrad`/`simplex`/...) | 0 | The C++ class is an abstract dispatch over (SeedGenerator × Builder); Julia uses multiple dispatch on FCN type, no inheritance needed. |
| `VariableMetricMinimizer.h` | (delete; not needed in Julia layout) | 0 | Just a (SeedGenerator, VariableMetricBuilder) bundle. |
| `VariableMetricBuilder.h/.cxx` | `JuMinuit._migrad_iterate!(workspace, fcn, seed, strategy, maxfcn, edmval)` (internal) | 0 | The 200-line MIGRAD loop in `VariableMetricBuilder.cxx:205–375`. Line-by-line port. |
| `MnSeedGenerator.h/.cxx` | `JuMinuit._seed_state(fcn, x0, errors, strategy)` | 0 | `MnSeedGenerator.cxx:42–101`. |
| `MnLineSearch.h/.cxx` | `JuMinuit._line_search!(workspace, fcn, ...)` | 0 | Parabolic interpolation; `MnLineSearch.cxx:46–313`. Cubic and Brent variants are conditional `#ifdef USE_OTHER_LS` — **do not port** (see Deferred §9). |
| `MnParabola.h/.cxx` + `MnParabolaFactory.h/.cxx` + `MnParabolaPoint.h` | inline helpers in `linesearch.jl` | 0 | Three tiny classes, fuse into local helpers. |
| `Numerical2PGradientCalculator.h/.cxx` | `JuMinuit._numerical_gradient!(workspace, fcn, p, prev_grad, strategy)` | 0 | The two-point central-difference algorithm at `Numerical2PGradientCalculator.cxx:63–230`. OpenMP block (line 112) → Phase 2.2. |
| `InitialGradientCalculator.h/.cxx` | `JuMinuit._initial_gradient!(workspace, p, fcn, trafo, strategy)` | 0 | Uses initial parameter step sizes (`MinuitParameter::Error`). |
| `HessianGradientCalculator.h/.cxx` | `JuMinuit._hessian_gradient!(...)` | 1 | Refined gradient inside `MnHesse`. |
| `AnalyticalGradientCalculator.h/.cxx` | `JuMinuit._analytical_gradient!(...)` (consumes user-provided `∇fcn`) | 1 / 2 | Phase 1 supports user-supplied gradient; Phase 2.1 wires up AD. |
| `DavidonErrorUpdator.h/.cxx` | `JuMinuit._davidon_update!(workspace, V, p1, g1, s0)` | 0 | `DavidonErrorUpdator.cxx:24–73`. Includes rank-1 branch when `delgam > gvg`. |
| `BFGSErrorUpdator.h/.cxx` | `JuMinuit._bfgs_update!(...)` | 2 | Alternative updator; `MnMigrad(BFGSType{})` path. Phase 2 (opt-in). |
| `FumiliBuilder/FumiliMinimizer/Fumili*` (10 files) | (deferred) | — | See Deferred §9. |
| `MinimumBuilder.h/.cxx` | folded into builder functions | 0 | The C++ abstract base only holds print level + tracer + storage level — represent in Julia as fields on the `Builder` callable struct. |
| `MinimumSeed.h` + `BasicMinimumSeed.h` | `MinimumSeed` (immutable struct) | 0 | No shared-ptr indirection needed. |
| `MinimumState.h` + `BasicMinimumState.h` | `MinimumState` (immutable struct, but builder mutates a current-state ref) | 0 | C++ uses `shared_ptr<BasicMinimumState>` (`MinimumState.h:66`); Julia uses a plain `struct` value. |
| `MinimumParameters.h` + `BasicMinimumParameters.h` | `MinimumParameters` (struct: `x::Vector{Float64}`, `dirin::Vector{Float64}`, `fval::Float64`) | 0 | |
| `MinimumError.h` + `BasicMinimumError.h` | `MinimumError` (struct: `inv_hessian::Symmetric{...}`, `dcovar::Float64`, status flags) | 0 | C++ status enums (`MnHesseFailed`, `MnMadePosDef`, etc.) → Julia `@enum CovStatus` or `Symbol`. |
| `FunctionGradient.h` + `BasicFunctionGradient.h` | `FunctionGradient` (struct: `grad::Vector`, `g2::Vector`, `gstep::Vector`, `analytical::Bool`) | 0 | |
| `FunctionMinimum.h` + `BasicFunctionMinimum.h` | `FunctionMinimum` (struct: seed, states, up, status flags) | 0 / 1 | P0 keeps only last state; P1 expands to full history if `storage_level==1`. |
| `MinimumErrorUpdator.h` | (abstract interface; Julia uses dispatch on updator type) | 0 | |
| `FCNBase.h` + `FCNGradientBase.h` + `FCNAdapter.h` + `FCNGradAdapter.h` | `CostFunction{F, T}` and `CostFunctionWithGradient{F, G, T}` (parametric structs) | 0 / 1 | The C++ inheritance hierarchy collapses to two concrete types; multiple dispatch handles the grad-vs-no-grad branching. |
| `MnFcn.h/.cxx` | call-counting wrapper inside the builder (`workspace.nfcn += 1` after each FCN call) | 0 | C++ wraps FCN to count calls; in Julia we increment in the workspace. |
| `MnUserFcn.h/.cxx` | folded into `MigradWorkspace.eval_fcn(x_internal)` which applies the internal→external transform | 1 | Phase 0 has no bounds, so this collapses to a direct FCN call. |
| `MnUserParameters.h/.cxx` | `Parameters` (struct holding `Vector{MinuitParameter}` + name map) | 1 | Phase 0 takes raw `Vector{Float64}`; P1 adds the named-parameter facade. |
| `MnUserParameterState.h/.cxx` | `ParameterState` (struct) | 1 | |
| `MnUserCovariance.h` | `Covariance` (struct over `Symmetric` or packed `Vector{Float64}`) | 1 | |
| `MnUserTransformation.h/.cxx` | `Transformation` (struct holding parameter metadata + cached external values) | 1 | Heart of the bounds path (`MnUserTransformation.cxx:99–141`). |
| `MinuitParameter.h` | `MinuitParameter` (struct) | 0 / 1 | Phase 0: just `value`, `error`, `fixed::Bool`. Phase 1: adds limits. |
| `MnGlobalCorrelationCoeff.h` | `GlobalCorrelation` (computed on demand) | 1 | |
| `SinParameterTransformation.cxx` + `SqrtUpParameterTransformation.cxx` + `SqrtLowParameterTransformation.cxx` | three `int2ext` / `ext2int` / `dint2ext` methods in `transform.jl` | 1 | Pure scalar math; 3 × ~50-line files. |
| `MnPosDef.h/.cxx` | `_make_posdef!(error, prec)` | 0 | Adds-to-diagonal trick to enforce positive definiteness; `MnPosDef.cxx:30–104`. |
| `MnMachinePrecision.h/.cxx` | `MachinePrecision` (struct: `eps`, `eps2`); default from `eps(Float64)` | 0 | |
| `MnStrategy.h/.cxx` | `Strategy` (struct); `Strategy(0)`, `Strategy(1)`, `Strategy(2)` constructors | 0 | Values from `MnStrategy.cxx:33–70`. |
| `VariableMetricEDMEstimator.h/.cxx` | `_estimate_edm(grad, error)` (free function) | 0 | One-liner: `0.5 * similarity(grad, inv_hessian)`. |
| `NegativeG2LineSearch.h/.cxx` | `_negative_g2_line_search!(...)` | 0 / 1 | Phase 0 may skip on the assumption initial g2 is positive on the chosen benchmarks; **flag the assumption in the seed code** and add it in Phase 1 with regression test. |
| `MPIProcess.h/.cxx` | (deferred) | — | See Deferred §9. |
| `LASymMatrix.h` + `LAVector.h` + `MnMatrix.h` | use `Symmetric{Float64,Matrix{Float64}}` + `Vector{Float64}` directly | 0 | Or packed buffers + custom `*spmv!` if benchmarking demands; decision in §3.5 day 4–7. |
| `ABObj.h` + `ABSum.h` + `ABProd.h` + `ABTypes.h` + `LaSum.h` + `LaProd.h` + `LaOuterProduct.h` + `LaInverse.h` + `LaProd.h` + `VectorOuterProduct.h` + `MatrixInverse.h` | (do not port; use `mul!`/`axpy!`/`syr!`) | — | The expression-template layer is a C++ workaround for what Julia handles natively. |
| `StackAllocator.h` | (do not port) | — | Workspace preallocation in `MigradWorkspace` is the Julia analog. |
| `MnRefCountedPointer.h` + `MnReferenceCounter.h` | (do not port) | — | Pre-`shared_ptr` reference counting; superseded by `shared_ptr` in C++ and irrelevant in Julia. |
| `mndaxpy.cxx` + `mndscal.cxx` + `mnddot.cxx` + `mndspmv.cxx` + `mndspr.cxx` + `mndasum.cxx` + `mnlsame.cxx` + `mnxerbla.cxx` | (do not port) | — | These are f2c-translated BLAS routines. Use `LinearAlgebra.BLAS.*` directly. |
| `mnvert.cxx` | `_sym_invert!(M)` calling `LAPACK.sptrf!`/`sptri!`, fallback to Gauss-Jordan for tiny n | 0 / 1 | Worth a unit test comparing both paths on a 4×4 case. |
| `mnteigen.cxx` + `LaEigenValues.cxx` | `_sym_eigvals(M)` via `LAPACK.spev!('N', 'U', ...)` or `eigvals(Symmetric(M))` | 0 / 1 | Needed inside `MnPosDef`. |
| `mntplot.cxx` + `MnPlot.h/.cxx` + `mnbins.cxx` | (use RecipesBase recipes instead) | 2 | C++ text-plot is for terminal use; Julia users have Plots/Makie. |
| `MnPrint.h/.cxx` + `MnPrintImpl.cxx` | use `Logging` stdlib + `@debug`/`@info` | 0 | C++ has a hand-rolled print system; Julia's logging is sufficient. |
| `MnTraceObject.h/.cxx` + `TMinuit2TraceObject.cxx` | optional callback in `migrad(... ; trace=callback)` | 1 / 2 | Trace each iteration; useful for debugging numerical-divergence cases. |
| `MnCovarianceSqueeze.cxx` | inline helper in `state.jl` | 1 | Strips fixed parameters from covariance matrix. |
| `MnEigen.h/.cxx` | `eigen(::Covariance)` method | 1 | Trivial. |
| `MinimizerOptions.cxx` (in `src/math/`) | `MinimizerOptions` keyword-args on `migrad` | 1 | |
| `ParametricFunction.h/.cxx` | (skip) | — | ROOT IFunction integration; not needed standalone. |
| `examples/simple/` | mirror as `examples/quad4.jl`, `examples/gauss_ll.jl` | 1 | Doc-driving examples. |
| `test/MnSim/*` + `test/MnTutorial/*` | mirror corresponding C++ tests as Julia tests; one Julia test = one C++ test | 0–1 | The Phase-0 corpus is `Quad1F`, `Quad4F`. |
| `Math/Minimizer.h`, `Math/IFunction*.h`, `Fit/ParameterSettings.h` | (skip) | — | ROOT framework glue; Julia uses native types. |

**Counting**: ~26K LOC C++ → expected ~5–7K LOC Julia (much of the C++ is
boilerplate around `shared_ptr` payloads, the expression-template layer,
and BLAS translations). Phase 0 alone should be < 1.5K LOC Julia + < 1K
test LOC.

---

## 8. Risk register

| # | Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|---|
| 1 | **DFP Hessian update numerical drift** vs C++. The rank-2/rank-1 branch in `DavidonErrorUpdator.cxx:62` is sensitive to summation order. Even bit-identical math through Julia BLAS vs reference BLAS will diverge over many iterations. | High | High — breaks the 1e-10 acceptance criterion on long fits. | (a) Compare iteration-by-iteration `inv_hessian` to a captured C++ trace on Rosenbrock-10 for ≤ 50 iter; flag the first iteration where any element diverges past 1e-12 and root-cause. (b) Document the divergence point in the test; accept that final answers agree to 1e-10 even if iteration counts drift by ±2. (c) Implement the matrix updates in exactly the same FLOP order as `DavidonErrorUpdator.cxx:58–69`. |
| 2 | **BLAS thread interaction** when a user wires multi-threaded likelihoods (Phase 2). Nested OpenMP-in-OpenBLAS deadlocks or wastes cycles. | Medium | Medium | Default `BLAS.set_num_threads(1)` when `Threads.nthreads() > 1`; document; provide a `with_blas_threads` helper. |
| 3 | **ABObj / expression-template substitute** under-performs. If our Julia replacement (broadcasting + `mul!`) misses fusion opportunities, the gradient-loop inner product chain (e.g. `v0 * dg` followed by outer products) might cost 2× extra. | Medium | High at Phase 0 gate | Profile early (Phase 0 day 8–12). If inner-loop allocations or extra passes appear, write fused micro-kernels for the 3–4 specific compound ops the MIGRAD loop uses (e.g. `step = -V*g`, `vg = V*dg`, `vUpd = αoouter(dx) − βouter(vg)`). |
| 4 | **Closure non-specialization**. If the user passes an `fcn::Function`-typed callable rather than allowing parametric specialization (`F`), every FCN call goes through dynamic dispatch. iminuit users may pass closures from notebooks. | Medium | High | Force specialization by holding the FCN as a parametric field in `CostFunction{F}` and by making every internal builder function generic on `F`. Use `Test.@inferred` and `JET.@report_call` in CI on a representative path. |
| 5 | **Numerical instabilities with bounded parameters**. The sin transform (`SinParameterTransformation.cxx:38`) clamps internal values to `[-π/2, π/2)` minus a margin. Float64 vs C++ may handle the boundary case slightly differently, drifting MIGRAD iteration counts on bounded fits. | Medium | Medium | Phase 1 scope; add a dedicated "stress" test fitting a Gaussian where the parameter starts at the bound. Use the same `prec.Eps2()` formula as C++ (`MnMachinePrecision.h:41`). |
| 6 | **Reference data generation cost**. Building the C++ benchmark/reference binaries is itself a multi-hour task (CMake + Minuit2 standalone). | Medium | Low | Commit JSON reference dumps under `test/reference_data/`; CI does not require the C++ build. Document the regen procedure under `tools/regen_reference.md`. Cap at 10 reference cases. |
| 7 | **MnPosDef eigenvalue path divergence**. When the Hessian goes non-pos-def, `MnPosDef.cxx:80` calls `eigenvalues` then adds to the diagonal. Different eigenvalue routines (LAPACK `spev` vs Julia's `LinearAlgebra`) may pick a different perturbation. | Low–Medium | Medium | Use `LAPACK.spev!` directly to match C++'s eigensolver choice. If still divergent, port the C++ Jacobi (`mnteigen.cxx`) verbatim for the pos-def branch only. |
| 8 | **API churn between Phases 0 → 1 → 3**. A Phase 0 user shouldn't see breaking changes. But the iminuit-style `Minuit(...)` constructor of Phase 3 is structurally different from the free-function `migrad(...)` of Phase 0. | Medium | Low–Medium | Phase 0 ships `JuMinuit.migrad` as the only public function. Phase 1 adds `Migrad`/`Hesse`/`Minos` types. Phase 3 adds `Minuit` as a new symbol. Nothing in Phase 0 is removed. Document this commitment in `CHANGELOG.md`. |

---

## 9. Deferred

Listed explicitly so future maintainers know these *are* known, not forgotten:

- **Fumili minimizer** (`FumiliBuilder`, `FumiliMinimizer`, `FumiliErrorUpdator`,
  `FumiliGradientCalculator`, `FumiliChi2FCN`, `FumiliMaximumLikelihoodFCN`,
  `FumiliStandardChi2FCN`, `FumiliStandardMaximumLikelihoodFCN`,
  `FumiliFCNBase`, `FumiliFCNAdapter`, `MnFumiliMinimize.h/.cxx`). Specialized
  for chi² / max-likelihood with Jacobian-style updates; useful for some HEP
  fits. Out of scope until there's a concrete user demand. ~3K LOC saved.
- **MPI support** (`MPIProcess.h/.cxx`, `MPI_SYNCH_PROC` guards in
  `Numerical2PGradientCalculator.cxx:102–214` and `MnHesse.cxx:240`). Replaced
  by `Distributed.jl` if/when needed in Phase 2+. Not for the v1.0 release.
- **BFGS Hessian updator** (`BFGSErrorUpdator.h/.cxx`,
  `VariableMetricMinimizer(BFGSType)`). Use Davidon (DFP, the Minuit default).
  BFGS as a Phase 2 toggle.
- **CombinedMinimizer / ScanMinimizer** (`CombinedMinimizer.h/.cxx`,
  `CombinedMinimumBuilder.cxx`, `ScanMinimizer.h`, `ScanBuilder.cxx`). The
  "Combined" path is essentially MIGRAD + SIMPLEX fallback; provide as a
  Julia composition (`migrad(...) || simplex(...)`) rather than a new
  minimizer type. Phase 2.
- **`MnLineSearch::CubicSearch` and `MnLineSearch::BrentSearch`** at
  `MnLineSearch.cxx:321–820`. These are `#ifdef USE_OTHER_LS` paths,
  disabled by default in the C++ build. The default parabolic search is
  what every Minuit2 user runs. Don't port unless explicitly demanded.
- **ROOT serialization compatibility** (`G__DICTIONARY` in
  `FunctionMinimum.h:18–20`, `LinkDef.h`, the entire `inc/Math` and
  `inc/Fit` headers). Not relevant outside a ROOT process.
- **`MnTinyMain` / FORTRAN-era utility code** (`MnTiny.h/.cxx`,
  `TMinuit2TraceObject.cxx`). Compatibility shims for the old Fortran
  Minuit; obsolete.
- **OpenMP at the gradient level** (`#pragma omp parallel`/`for` blocks in
  `Numerical2PGradientCalculator.cxx`). Replaced by `Threads.@threads` in
  Phase 2.2.
- **`StackAllocator`** — moot in Julia (see Risk #3 mitigation).
- **`ABObj` and the entire expression-template layer** — replaced wholesale by
  Julia's broadcasting + LAPACK calls (see Module mapping §7).
- **Hand-rolled BLAS** (`mndaxpy`, `mndscal`, `mnddot`, `mndspmv`, `mndspr`,
  `mndasum`, `mnlsame`, `mnxerbla`) — use `LinearAlgebra.BLAS.*` directly.
- **`MnPlot` text plotting** (`MnPlot.cxx`, `mntplot.cxx`, `mnbins.cxx`) —
  Julia users get RecipesBase recipes in Phase 2.3 instead.
- **ParametricFunction integration** (`ParametricFunction.h/.cxx`) — ROOT
  function-object integration; users in Julia pass plain callables.

---

## 10. Open questions for the user

A handful of decisions need the user's judgment before code is written.

1. **Bounded-parameter representation**. iminuit and the C++ code expose
   *external* (bounded) values to the user and *internal* (unbounded, via
   sin/sqrt transform) to the optimizer. Internally we must mirror this for
   numerical equivalence. **Question**: should the public Julia API return
   covariance / errors in *external* coordinates only (mirroring iminuit /
   IMinuit.jl), or also expose the internal representation (mirroring the
   C++ `MnUserParameterState::IntCovariance()`) for advanced users?

2. **Bounds transform deviation**. The C++ uses sin (double-bounded) and sqrt
   (one-sided) transforms — these are *Minuit's* historical choices, not
   universally optimal. iminuit users in modern HEP have occasionally
   complained about the sin-transform pathology near the bounds (gradient → 0
   as the parameter approaches the limit). **Question**: do we strictly mirror
   the C++ transforms for numerical equivalence, or offer optional tanh /
   logistic alternatives behind a flag? (Strict mirror recommended for v1.0.)

3. **Default linear-algebra storage**. The C++ uses lower-triangular packed
   storage for symmetric matrices (`LASymMatrix`) for memory locality. Julia
   can do the same with `Vector{Float64}` + custom BLAS, *or* use
   `Symmetric{Float64,Matrix{Float64}}` with the dense matrix being 2× the
   memory but better cache behavior for small n. **Question**: any preference,
   or decide by benchmark in Phase 0 day 4–7?

4. **AD backend choice for Phase 2.1**. ForwardDiff.jl is mature and
   composes with all Float64 code; Enzyme.jl is faster on hot inner loops but
   has rougher edges. **Question**: pick ForwardDiff as primary and Enzyme as
   opt-in, or the reverse? (ForwardDiff-primary recommended for the v1.0 user
   base.)

5. **Threads as default**. iminuit and the C++ default to single-threaded
   (OpenMP is opt-in via CMake flag). Julia code increasingly assumes
   `Threads.nthreads() > 1` is normal. **Question**: should the Phase 2.2
   parallel gradient default to ON when `nthreads() > 1`, or stay opt-in?
   (Opt-in recommended for predictability and benchmarking parity.)

6. **API parity vs. idiomatic Julia at Phase 3**. iminuit's API uses Python
   property access (`m.values`, `m.errors`) and a builder-style fluent
   interface. **Question**: how literally should JuMinuit mirror this? Options:
   (a) `m.values` via `Base.getproperty` overloading (faithful but un-idiomatic),
   (b) `values(m)` / `errors(m)` accessors (idiomatic but breaks copy-paste),
   (c) both. (Recommend both; deprecate (a) in v2.0.)

7. **Status enums in `MinimumError`**. C++ uses tag types
   (`MnHesseFailed`, `MnMadePosDef`, `MnInvertFailed`, `MnNotPosDef`) for
   error-matrix status. Julia equivalents: `@enum CovStatus`, `Symbol`,
   or `Val{:hesse_failed}`. **Question**: pick one; recommend `@enum` for
   inferability and printability.

8. **Reference data freshness**. The C++ reference data we commit will need
   updating if we ever bump to a newer Minuit2 release. **Question**: pin
   to the current GooFit/Minuit2 commit as the canonical reference for v1.0,
   document the regen procedure, and address upgrades in v2.0?

---

## 11. Critical files for Phase 0 implementation

Three Julia files and one C++ file together define the Phase 0 surface and the
oracle against which everything is measured.

- `/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/src/migrad.jl` — the MIGRAD loop
- `/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/src/davidon.jl` — the DFP update (correctness-critical)
- `/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/src/gradient.jl` — the numerical gradient (performance-critical)
- `/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/src/linesearch.jl` — the parabolic line search
- `/Users/fkg/Coding/Agents/ResearchWork/JuMinuit/reference/Minuit2_cpp/src/VariableMetricBuilder.cxx` — line-by-line oracle for the Julia loop
```

Now I'll output this as the final response. Given the length constraint and the user's request to "return only the ROADMAP.md markdown text (no commentary around it)", I'll deliver the markdown directly. I also need to follow the system prompt's final section requirement about listing "Critical Files for Implementation" — but that section IS in the ROADMAP itself as §11. The system prompt says to end the response with that section listing 3-5 critical files, and the ROADMAP has that. So the ROADMAP body satisfies both the user request and the system prompt's structural requirement simultaneously.
