# Benchmark results — JuMinuit.jl vs IMinuit.jl on real physics fits

This file captures the most recent comparison runs of the two example
fits across all available execution schemes. Re-run with:

```bash
julia -t 8 --project=scripts BenchmarkExamples/X3872_dip/bench_full.jl
julia -t 8 --project=scripts BenchmarkExamples/IAM_2Pformfactor/bench_full.jl
```

Each script reports stage-by-stage median wall-time (3 rounds + warmup)
and performs cross-checks on every stage (minimum, MINOS errors,
mncontour centroid). The summary tables below reflect the latest run on
macOS / Julia 1.12 / `julia -t 8` / `BLAS.set_num_threads(1)`. See the
commit history for older runs.

## Scheme legend

| label       | description                                              |
|-------------|----------------------------------------------------------|
| `jm_num`    | JuMinuit numerical gradient, sequential                  |
| `jm_ad`     | JuMinuit AD (ForwardDiff) — package extension            |
| `jm_th_num` | JuMinuit threaded numerical (Phase G)                    |
| `jm_th_ad`  | JuMinuit threaded AD                                     |
| `iminuit`   | Python `iminuit` via PyCall (IMinuit.jl `v0.2.1`)        |

## X(3872) dip fit — 3 params, FCN ~ 38 μs/call, 4 data points

Published analysis: V. Baru, F.-K. Guo, C. Hanhart, A. Nefediev,
"How does the X(3872) show up in e⁺e⁻ collisions: Dip versus peak",
*Phys. Rev. D* **109** (2024) 11, L111501,
[arXiv:2404.12003](https://arxiv.org/abs/2404.12003),
[INSPIRE 2778938](https://inspirehep.net/literature/2778938).

| scheme       | migrad+hesse | minos (3 params) | mncontour (20 pts) |
|--------------|-------------:|-----------------:|-------------------:|
| `jm_ad`      |  **3.3 ms**  |  **75 ms**       |  **35 ms**         |
| `jm_th_ad`   |  3.5 ms      |  74 ms           |  35 ms             |
| `jm_num`     |  5.3 ms      | 129 ms           |  90 ms             |
| `jm_th_num`  |  5.9 ms      | 138 ms           |  97 ms             |
| `iminuit`    |  7.2 ms      | 155 ms           |  50 ms             |

(latest run on commit-tree post gap closure; within ±5 % of the
pre-closure baseline — no regression from the 9 merged PRs M1–P5).

**Headlines**

- JuMinuit AD is **2× faster than JuMinuit numerical** on migrad, **1.7×** on minos.
- JuMinuit numerical is **30 % faster than iminuit** on migrad, comparable on minos.
- JuMinuit AD reaches **2.3× iminuit speed** on migrad and **2.1×** on minos.
- 3-dim problem is too small for threading to help (`jm_th_*` ≈ sequential).

**MNCONTOUR caveat for X(3872)**

The X(3872) fit overfits 3 parameters on 4 points (χ²_min = 0.017), so
the 1σ region collapses to near machine precision in some directions.
Both libraries terminate early with "MnContours unable to find first
two points" on every parameter pair tested; the wall-times above are
the time spent until early termination, not the time of a successful
contour generation.

**Open issue (X3872)**

JuMinuit's MINOS returns `(0, 0)` for par[2] (`r`) at a minimum where
iminuit successfully returns `(-0.00214, +0.00431)`. Both backends
converge to the same x, fval (Δx ≈ 9·10⁻⁶, Δfval ≈ 9·10⁻⁸), so this
is a JuMinuit edge-case in `function_cross` for tight wells, not a
fit-quality artifact. Tracked.

## IAM 2π form-factor — 9 LECs, FCN ~ 10 ms/call, 85 data points

| scheme       | migrad+hesse | minos (par 1)        | mncontour (8 pts) |
|--------------|-------------:|---------------------:|------------------:|
| `jm_num`     | **5.42 s**   | **16.2 s**           | **27.3 s**        |
| `iminuit`    | 18.74 s      | REFUSED (invalid fmin) | REFUSED         |
| `jm_ad`      | FAILED       | —                    | —                 |
| `jm_th_*`    | SKIPPED (Phase H rejects) | —       | —                 |

(latest run on commit-tree post gap closure; within ±2 % of the
pre-closure baseline — no regression from the 9 merged PRs M1–P5).

**Headlines**

- JuMinuit MIGRAD is **3.5× faster than iminuit** on the 9-LEC fit
  (5.4 s vs 19.0 s) — **but lands at a SHALLOWER minimum** (fval=613.5
  vs fval=409.9). The no-improvement early-termination check in
  `src/migrad.jl` triggers too aggressively on this landscape. Both
  fits report `is_valid=false` (above-max-edm). Tracked as follow-up.
- iminuit hard-refuses MINOS / MNCONTOUR on an invalid fmin (Python
  raises `RuntimeError("Function minimum is not valid")`). JuMinuit
  runs both to completion on the same invalid fmin — MINOS returns
  `(0, 0)` for this tight well, MNCONTOUR returns an empty point set.
  Neither behavior is "correct"; both libraries struggle. The bench
  wraps the iminuit calls in `try/catch` so the script still completes
  end-to-end.
- **Phase H pre-flight catches IAM thread-unsafety in milliseconds**:
  `is_thread_safe(chi2_iam, paras0)` returns `false` because
  `St4_00!` writes a module-level `const c_00_4` buffer. All `jm_th_*`
  schemes are refused before any migrad work happens — this is the
  silent-wrong-answer fix from commit `96513d7` demonstrated on a
  real physics fit.
- AD path FAILS: IAM's `src/amplitudes.jl` etc. carry `Float64`
  annotations that block ForwardDiff `Dual` propagation. Genuine
  limitation of the IAM source, not a JuMinuit issue.

## Methodology

- **Wall-time**: 3 rounds (X3872: 5) + 1 warmup, take the median.
  `GC.gc()` between rounds. `sleep(0.2–0.5 s)` between rounds.
- **Cross-checks**: every stage compares all paths' results against
  `jm_num` (the most-conservative reference). Mismatches are flagged
  but do not abort the bench.
- **Phase H**: when `Threads.nthreads() > 1`, the bench probes
  `is_thread_safe(cf, x0)` before launching threaded schemes.
  Racey FCNs are refused upfront.
- **FCN cost**: measured with `@benchmark` (`BenchmarkTools.jl`),
  reported in the per-script header.

## Closed follow-up work

Both originally-flagged follow-ups were closed by a single PR:

1. ✅ **MINOS early termination on tight wells** — X(3872) par[2] now
   returns `(-0.00439, +0.00439)` instead of `(0, 0)`; IAM par1 returns
   `(-0.000173, +0.000173)` instead of `(0, 0)`. Resolved by PR #6
   (commit `a1fa015`): loosening the `_migrad_loop` seed-acceptance
   gate to match C++ `BasicMinimumSeed::IsValid()` semantics (the C++
   check only looks at the seed's own `fValid` flag, not state
   validity), so a `MnHesseFailed`-status seed with structurally
   valid params + gradient + diagonal V is no longer rejected.

2. ✅ **IAM 9-LEC early-termination divergence** — at Strategy(2), IAM
   MIGRAD now reaches **fval = 401.45** (was stuck at 613.49),
   matching iminuit Strategy(0)'s **400.23** (same basin). At
   Strategy(2), X(3872) drops from 1.30 to **0.017** (matches the
   published [arXiv:2404.12003](https://arxiv.org/abs/2404.12003)
   global minimum). Resolved by the same PR #6: dropping a C++
   Minuit2 bug in `_hesse_diagonal_failure` — the second `eps2` clamp
   on `1/g2` was inverting the truth (mapping `1/g2 = 1e-10`, i.e.
   "very well determined", to `1.0`, i.e. "poorly determined"),
   producing `V = I` whenever any parameter was FCN-flat. Verified
   that iminuit Strategy(2) hits the same trap with the C++ formula,
   so this is a real upstream bug that JuMinuit now fixes.

Note: the default Strategy(1) IAM behavior is unchanged from the
pre-PR-#6 numbers above (the fix is gated on the
seed-time-MnHesse-bootstrap path that only Strategy(2) takes today —
the `_migrad_loop` strategy-1 cold path doesn't go through MnHesse at
the seed). The headline IAM 9-LEC fitter should now set
`m.strategy = Strategy(2)` for the deeper minimum.
