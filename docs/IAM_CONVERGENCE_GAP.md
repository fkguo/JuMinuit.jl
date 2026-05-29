# IAM cold-start convergence gap: JuMinuit vs iminuit (613 → 330)

**Date**: 2026-05-29
**Branch**: `feat/iam-convergence-gap`
**FCN**: IAM 2π form-factor 9-LEC fit (`BenchmarkExamples/IAM_2Pformfactor`)
**Seed**: `paras0 = [lecr0…, 1e-4]`, `error = fill(1e-6, 9)`, default tol.

## Symptom (as reported)

`bench_full.jl` `build_jm_num()` does `Minuit(chi2_iam, paras0; error=errs0)`
→ `migrad!(m)` → `hesse(m)` and prints:

```
jm_num   migrad+hesse: fval=613.485  is_valid=false
iminuit  migrad+hesse: fval=409.885  is_valid=false
```

JuMinuit lands in a shallower basin (613) than iminuit (409). The robust
retry layer (`migrad!` default `iterate=5, use_simplex=true`) runs but its
fixed-point detector confirms it cycles at 613.

## TL;DR — root cause and fix

**The gap is a high-level default-strategy mismatch, not a core MIGRAD defect.**

- iminuit's `Minuit` class defaults to `strategy = 1`; C++ Minuit2's
  `MnStrategy()` default is also level 1 (`SetMediumStrategy`).
- JuMinuit's high-level `Minuit(fcn, x0)` constructor defaulted to
  `Strategy(0)` — a Phase-0 holdover (`strategy.jl` docstring still read
  *"Phase 0 default and only supported level"*). The bench used the default,
  so **JuMinuit ran Strategy 0 while iminuit ran Strategy 1**. They were never
  running the same algorithm.

At **matched strategy**, JuMinuit's MIGRAD is *not* worse — at Strategy 1
(iminuit's own default) JuMinuit's single MIGRAD reaches **330.75**, deeper
than iminuit's retry result of 409.89.

**Fix**: change the high-level `Minuit(...)` constructor default from
`Strategy(0)` to **`Strategy(1)`, uniformly** — for both numerical and
analytical/AD (`grad=`) FCNs. iminuit's `Minuit` class applies `strategy=1`
regardless of whether a gradient is supplied, so JuMinuit must too; a default
that differed by construction path would be the *same class* of silent
mismatch that caused this gap. Two points:

- **AD path extended to all strategy levels** — the AD `seed_state`
  (`ad_gradient.jl`) previously threw `Strategy(0) only` (a Phase-0 holdover).
  It now mirrors the numerical `seed_state`: the diagonal-from-g2 seed is
  strategy-independent (the AD gradient is exact), Strategy 2 adds the
  seed-time MnHesse bootstrap, and Strategy ≥ 1's inner-HESSE refinement runs
  in `_migrad_loop` — all via `hesse(::AbstractCostFunction)`, which
  finite-differences `cf.f` and already accepted `CostFunctionWithGradient`.
  So `Minuit(fcn, x0; grad=g); migrad!(m)` now runs at S=1 instead of
  throwing.
- **Low-level entry points unchanged** — `migrad(cf, …)` / `seed` /
  `function_cross` / `minos` / `contours` keep their own `Strategy(0)`
  defaults (pinned to the C++ oracle reference data; `test_cpp_oracle.jl`
  asserts `strategy_level == 0`).

Implemented as `strategy = Strategy(1)` in both `Minuit(...)` constructor
methods (`src/minuit.jl`), plus the AD-seed extension in `src/ad_gradient.jl`.

Result after fix (same bench call, default settings):

```
jm_num   migrad!(m): fval=330.753   (was 613.485)   ← now BEATS iminuit's 409.885
```

## The data

### 1. Per-strategy single-shot (`iterate=1`, true apples-to-apples)

Both libraries, IAM `paras0` cold start, retry **disabled** on both sides
(iminuit via `m.migrad(iterate=1)`, JuMinuit via `migrad!(m; iterate=1)`):

| Strategy | JuMinuit single-shot | iminuit single-shot |
|----------|---------------------:|--------------------:|
| **S=0**  | 613.49               | 476.15              |
| **S=1**  | **330.75**           | 614.95              |
| **S=2**  | 1268.65 (stuck)      | 1268.65 (stuck)     |

### 2. Per-strategy with each library's native retry

iminuit retry = re-run the *same* `MnMigrad` (same strategy) from the last
point up to 5×. JuMinuit retry = growing Simplex hop + `Strategy(2)` MIGRAD,
up to 5×, with fixed-point/saturation stops.

| Strategy (pass 1) | JuMinuit retry | iminuit retry |
|-------------------|---------------:|--------------:|
| **S=0**           | 613.49         | 400.23        |
| **S=1**           | **330.75**     | 409.89 (default) |
| **S=2**           | 1268.65        | 1268.65       |

**Defaults**: JuMinuit was `S=0 → 613`; iminuit is `S=1 → 409`. After the fix
JuMinuit is `S=1 → 330`.

### 3. Aligned per-iteration trace at matched strategy (S=0)

Captured with `print_level=3`. Both seeds are identical:

```
                       JuMinuit S=0            iminuit S=0
  seed: FCN            1268.645892             1268.645892
  seed: Edm            1026.833019             1026.833283
  initial grad x0      -5.10e7  (g2 5.61e12)   -5.10e7  (g2 5.61e12)
  initial grad x8(=p9)  0       (g2 0)          0       (g2 0)   ← FCN ignores pars[9]
```

(`chi2_iam` reads only `pars[1:8]`; the 9th LEC is a flat direction. Both
libraries detect Negative-G2 on the seed and freeze it at 1e-4 — no
divergence here.)

First 10 DFP iterations — **identical to ~6 significant figures**:

| iter | JuMinuit S=0 fval | iminuit S=0 fval |
|-----:|------------------:|-----------------:|
| 0 | 1268.645892 | 1268.645892 |
| 1 | 987.2393689 | 987.2393917 |
| 2 | 978.3166239 | 978.3166157 |
| 3 | 962.6234493 | 962.6234307 |
| 4 | 930.1469850 | 930.1468527 |
| 5 | 912.4470956 | 912.4470872 |
| 6 | 909.5137347 | 909.5137144 |
| 7 | 907.3719752 | 907.3720844 |
| 8 | 896.3723745 | 896.3712664 |
| 9 | 887.5449898 | 887.5459662 |

The DFP update, EDM estimator (`edm·(1+3·dcovar)`), line search (accepted
step lengths α = 0.37, 0.15, 1.39, 9.95, 1.97, …), and Negative-G2 seed
handling are **byte-for-byte equivalent** to C++/iminuit for the first ten
steps. The divergence appears only later (≈ iter 20+) as the coarse
Strategy-0 2-cycle gradient accumulates noise into the two DFP trajectories.

**Localized divergence**: JuMinuit S=0 bails at **iter 24, fval = 613.49**
with *"matrix not pos.def; MnPosDef applied"* — the trial step gives
`gdel = step·g > 0` even after `MnPosDef`, so `_migrad_loop` breaks
(migrad.jl step 2, mirroring C++ `VariableMetricBuilder.cxx`). iminuit's S=0
first pass instead reaches 476, then its re-seed retry walks to 400. At
**Strategy 1**, JuMinuit's `dcovar`-triggered inner-HESSE refinement re-seeds
the curvature mid-run (at iter 49, `dcovar` climbs back to 0.12 → inner
HESSE fires) and the DFP loop descends all the way to **330.75** before
terminating — the path iminuit needs five retry passes to approximate.

## Hypothesis classification (the four candidates)

| Hypothesis | Verdict | Evidence |
|------------|---------|----------|
| **#4 Strategy-default mismatch** | ✅ **ROOT CAUSE** | JuMinuit default `Strategy(0)`, iminuit/C++ default level 1. Bench ran S=0 (→613) vs iminuit S=1 (→409). At matched S=1, JuMinuit (330) beats iminuit (409). |
| **#3 Seed inverse-Hessian scale** | ❌ ruled out | Seed Edm matches (1026.833019 vs 1026.833283); initial g, g2 identical per parameter; first DFP step lands at the same fval (987.239) on both sides. |
| **#1 Line-search timidity** | ❌ ruled out | Accepted step lengths are healthy/large (α up to 67.8), and the per-iteration fval matches iminuit step-for-step for the first 10 iters — the line search takes the *same* steps as C++. |
| **#2 No-improvement early exit** | ❌ ruled out as primary | JuMinuit ran 24 (S=0) / 49 (S=1) iterations before terminating — it did not exit on the no-improvement test (`|Δf| ≤ |f|·eps`). The S=0 termination at iter 24 is a `gdel>0`-after-MnPosDef pos-def bail, a downstream consequence of coarse-S0-gradient noise, not a too-tight threshold. At S=1 it does not bail there. |

## Why the fix is correct (not an IAM special-case)

1. **Drop-in fidelity** — the project's stated goal (commit `a884742`,
   "IMinuit.jl drop-in"). iminuit's `Minuit(fcn, x0)` defaults `strategy=1`;
   JuMinuit's high-level constructor must too, so a bare `migrad!(m)` matches
   `m.migrad()`. The `Strategy(0)` default was explicitly a Phase-0 limitation
   that was never updated when `hesse.jl` shipped (it enabled S≥1).
2. **It is a global default change**, not an IAM branch. Every high-level fit
   now runs the more thorough Strategy-1 path by default, exactly as iminuit.
3. **It satisfies the success criterion** — single-shot MIGRAD at the new
   default reaches 330.75 ≤ 410 on IAM.
4. **Low-level oracle parity preserved** — `migrad(cf, …)` etc. keep
   `Strategy(0)`; `test_cpp_oracle.jl` (which uses the low-level API and pins
   `strategy_level == 0`) is unaffected.

## Secondary findings (documented, not fixed here)

These are real and worth recording, but they live at a *non-default* strategy
after the fix, so they do not affect the headline bench. They are out of scope
for this P0 (changing retry semantics risks the X(3872) multistart feature
PR #8/#12 added).

- **At Strategy 0, JuMinuit's single MIGRAD (613) underperforms iminuit's
  (476).** Cause: the Strategy-0 2-cycle numerical gradient is noisy; the DFP
  trajectory diverges from iminuit's after ~20 iters and bails on a
  non-pos-def trial step. This is inherent to "fast/loose" S=0 and is why
  neither library validates at S=0 single-shot.

- **JuMinuit's retry differs from iminuit's, and is less effective at S=0.**
  iminuit's `Minuit.migrad` (verified by reading its source) retries by
  re-running the *same* `MnMigrad` object — **same strategy, no Simplex hop,
  no strategy bump** — from the last point, which freshly re-seeds the
  diagonal curvature each pass (476 → 400 at S=0). JuMinuit's retry instead
  bumps numerical FCNs to `Strategy(2)` (migrad! / minuit.jl `retry_strategy`)
  and does a growing Simplex hop. The in-code comment there claims
  *"Strategy(2) … the iminuit default"* — **this is inaccurate**; iminuit does
  not bump strategy on retry. The bump pushes a cold seed into the S=2
  pathology (see below), so JuMinuit's S=0 retry cannot match iminuit's 400.
  After the default fix this is moot for the default path (pass 1 is S=1 →
  330; the doomed S=2 retry is caught by fixed-point detection at
  `n_passes=2`), but the comment should be corrected and the bump
  reconsidered in a follow-up.

- **Strategy 2 from a cold seed is pathological for *both* libraries** (both
  stuck at 1268.65 — exact parity). At `paras0` the gradient is ~1e6 and the
  MnHesse-fail fallback yields `V ≈ I` (the C++ second clamp restored in
  PR #10), so the first Newton step `−V·g` has magnitude ~1e6 and the line
  search cannot make progress. This confirms the
  `docs/DAVIDON_CXX_AUDIT.md` "S=2 cold-seed pathology" note and is *not* a
  JuMinuit bug.

## Verification

- IAM `paras0`, default `migrad!(m)`: **613.49 → 330.75** (beats iminuit 409.89),
  `n_passes=2`. Guarded by
  `BenchmarkExamples/IAM_2Pformfactor/test_convergence_gap.jl`.
- IAM `paras0`, single-shot `migrad!(m; iterate=1)` at the default: 330.75 ≤ 410 ✓.
- Single-shot S=0/S=2 and all iminuit numbers **unchanged** (fix touches only
  the high-level default, not the algorithm).
- New unit regression `test/test_minuit.jl::"Default strategy = 1 (iminuit
  Minuit-class parity)"` asserts: numerical default → `Strategy(1)`, AD
  (`grad=`) default **also → `Strategy(1)`** and a default AD fit runs
  end-to-end at S=1 (plus an AD-at-S=2 fit), explicit strategy respected.
- `test/test_minuit_retry.jl::"Retry actually triggers …"` was pinned to
  `Strategy(0)` (it exercises the *retry* branches via a pass-1 stall that the
  coarse S=0 gradient reliably produces; the new S=1 default descends too far
  for the retry to enter, which would leave those branches uncovered). The
  AD-retry testset's stale "seed rejects level != 0" comment was corrected.
  No other test changed.
- Full `Pkg.test()` passes (see PR). The `_hesse_diagonal_failure` clamp +
  do-while loop (PR #10, which gets the IAM x_jm warm start to 322 via the S=2
  retry path) are untouched: that path goes through the low-level `migrad` and
  the retry S=2 bump, neither of which this change modifies. `test_cpp_oracle.jl`
  (low-level `migrad`, `strategy_level == 0`) is likewise unaffected.

## Relationship to `docs/DAVIDON_CXX_AUDIT.md`

That audit studied the **x_jm warm start** (S=2 basin walk to 322) and recorded
in its post-resolution table:

> | paras0 S=1 cold-shot | … | 330.75 (x_jm 8D basin + x[9]=1e-4) | 409.89 (x_im basin) |

i.e. it already *measured* that JuMinuit at S=1 reaches 330.75 and iminuit at
S=1 reaches 409.89 — but it did not connect that to the constructor's
`Strategy(0)` default, because its focus was the warm-start regime. This
document closes that loop: the cold-start bench gap is the default-strategy
mismatch, and matching iminuit's default both closes and reverses it.
