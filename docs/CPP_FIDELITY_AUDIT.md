# C++ Minuit2 ↔ JuMinuit line-by-line fidelity audit

**Date**: 2026-05-29 · **Branch**: `feat/iam-convergence-gap`
**Reference**: `reference/Minuit2_cpp/{src,inc}/*` (GooFit/Minuit2 v6.24.0)
**Scope**: deep, branch-by-branch comparison of individual ported algorithms —
*not* the component-level coverage map (that lives in `DEFERRED.md` /
`GAP_AUDIT.md`). Each section maps every C++ branch/exit path to its JuMinuit
counterpart and classifies it: ✓ faithful · documented-divergence · minor ·
missing.

Audited so far:

1. [MnHesse](#1-mnhesse) — `MnHesse.cxx:93-316` ↔ `src/hesse.jl`
2. [VariableMetricBuilder / MIGRAD](#2-variablemetricbuilder--migrad) —
   `VariableMetricBuilder.cxx` ↔ `src/migrad.jl:_migrad_loop`
3. [MnMinos](#3-mnminos) — `MnMinos.cxx` (+ `MnFunctionCross.cxx`) ↔
   `src/minos.jl` / `src/function_cross.jl`

---

## 1. MnHesse

`MnHesse.cxx:93-316` (the `operator()(MnFcn, MinimumState, MnUserTransformation,
maxcalls)` "real Hessian calculation"). Lines 318-414 are dead commented-out
code, ignored.

### Branch map

| C++ (MnHesse.cxx) | JuMinuit (hesse.jl) | Verdict |
|---|---|---|
| `amin=mfcn()`, `aimsag=√eps2·(\|amin\|+Up)`, `maxcalls=200+100n+5n²` (102–109) | 96–97, 91–93 | ✓ |
| init `g2/gst/grd/dirin=gst/yy` (112–116) | 108–112 | ✓ |
| analytical-gradient g2/step recompute (120–126) | 166–180 | ✓ (2 documented nuances) |
| diagonal `dmin=8·eps2·(\|xtf\|+eps2)`, `d=\|gst\|` (136–139) | 192–194 | ✓ |
| 5× multiplier loop, `sag≠0→break` (147–169) | 205–221 | ✓ except limits branch |
| L26 sag-zero → diagonal fallback `MnHesseFailed` (171–183) | 223–226 | ✓ |
| L30 `g2=2·sag/d²`, `grd`, `d=√(2·aimsag/\|g2\|)` (185–197) | 228–238 | ✓ except `d=min(0.5,d)` limits clamp |
| convergence `Tolerstp`/`TolerG2`, `d∈[0.1,10]·dlast` (203–208) | 241–256 | ✓ (defensive `g2≠0` guard, same result) |
| `vhmat(i,i)=g2(i)` (210) | 259 | ✓ |
| maxcalls-exhausted → diagonal fallback (211–223) | 269–275 | ✓ |
| Strategy>0 HGC gradient refine (228–235) | 290–303 | ✓ |
| off-diagonal `(fs1+amin−yy_i−yy_j)/(dirin_i·dirin_j)` (239–272) | 307–329 | ✓ (simple `i<j` = C++'s own old form) |
| `MnPosDef` on H (278) | 342 | ✓ (passes H not V — matches C++) |
| `Invert`; fail → diagonal fallback `MnInvertFailed` (283–296) | 348–355 | ✓ |
| `IsMadePosDef` → `MnMadePosDef` state (302–306) | 359–364 | ✓ |
| accurate → `dcovar=0` state (309–315) | 358–375 | ✓ |
| double-clamp `g2<eps2?1:1/g2; <eps2?1` ×3 fallbacks (177–180/216–219/289–292) | `_hesse_diagonal_failure` 462–463 | ✓ (abs-variant, identical result) |
| MPI off-diagonal partitioning (240–271) | — | intentionally not ported (MPI deferred) |

### Findings

- **MISSING (documented, narrow): bounded-parameter step clamping.**
  `has_limits = false` is hardcoded (hesse.jl:187), so two C++ branches never
  fire: the multiplier-loop `if HasLimits && d>0.5 → d=0.51`/fail (160–167) and
  the L30 `if HasLimits → d=min(0.5,d)` (194–195). HESSE runs in internal
  (arcsin) coordinates, where C++ clamps the probe step `d≤0.5` for
  externally-bounded params (near a bound the transform is steep; an unclamped
  `d` → wild external excursion → wrong 2nd-derivative). **Severity NICE-TO-HAVE**:
  only binds for poorly-determined / near-bound params or flat directions that
  trip the ×10 multiplier; well-determined bounded params away from bounds keep
  `d` small so it never fires. Documented as a Phase-1 first-cut deferral
  (hesse.jl:184–187; `DEFERRED.md` "bounds integration is the follow-up").
  Unbounded fits are fully faithful.

- **Documented faithful-but-different (not gaps):**
  - *Analytical-gradient gate* (hesse.jl:150–165): gated on `cf isa
    CostFunctionWithGradient` vs C++'s `IsAnalytical()` flag → a repeat `hesse`
    call re-refreshes (idempotent, extra FCN calls only, not a correctness bug).
  - *Analytical seed semantics* (132–148): recompute seeds from stale
    `state.gradient` vs C++'s fresh per-parameter user errors
    (`InitialGradientCalculator`); converges identically for smooth FCNs, can
    differ for pathological ones (GAP_AUDIT P2 follow-up).
  - *`abs()` in the double-clamp*: same result as C++'s raw comparisons
    (negative g2 → 1.0 both ways).
  - *Off-diagonal loop*: simple nested `i<j` vs C++'s MPI-flattened index
    arithmetic — mathematically identical (it *is* C++'s own non-MPI form,
    lines 400–410 of the commented block).

**Verdict: faithful port.** Every branch, exit path, formula, tolerance, and
the load-bearing double-clamp are correct. The only real omission is the
bounded-parameter step clamp — documented, narrow, unbounded-fits-unaffected.

---

## 2. VariableMetricBuilder / MIGRAD

`VariableMetricBuilder.cxx` ↔ `src/migrad.jl:_migrad_loop`. C++ splits MIGRAD
into an **outer** `Minimum` (54–203: edmval scaling, validity gates, the
do-while calling the inner loop + Strategy≥1 HESSE refinement) and an **inner**
`Minimum` (205–375: the DFP iteration). JuMinuit inlines both into one
`_migrad_loop` (outer `while iterate` wrapping inner `while true`) — same
control flow.

### Inner DFP loop (C++ 205–375 ↔ migrad.jl 690–878)

| C++ branch | JuMinuit | Verdict |
|---|---|---|
| `edm *= (1+3·Dcovar)` (229) | 586, 844 | ✓ |
| `step = −V·g` (241) | 724 `sym_mul!` | ✓ |
| zero-grad `⟨g,g⟩≤0 → break` (247–250) | 727–729 | ✓ |
| `gdel = step·g` (252) | 731 | ✓ |
| `gdel>0` → MnPosDef → recompute → still>0 → exit (254–273) | 734–748 | ✓ |
| line search (275) | 752 `line_search` | ✓ |
| no-improvement `\|pp.Y−Fval\|≤\|Fval\|·Eps → break` (278–291) | 762–767 | ✓ (≤eps·\|fval\| micro-diff) |
| accept `p = x + pp.X·step` (296) | 778 | ✓ |
| new grad `g = gc(p, s0.grad)` (298) | 785 | ✓ |
| `edm = Estimate(g, s0.Error())` — OLD error (300) | 792 | ✓ |
| `isnan(edm) → break` (302–306) | 794–796 | ✓ |
| `edm<0` → MnPosDef → recompute → still<0 → exit (308–321) | 799–806 | ✓ |
| Davidon `Update(s0,p,g)` (322) | 834–840 | ✓ |
| `while edm>edmval && nfcn<maxfcn` (341) | 878 | ✓ |

### Outer loop + finalization (C++ 54–203 ↔ migrad.jl 530–973)

| C++ branch | JuMinuit | Verdict |
|---|---|---|
| `edmval *= 0.002` (66) + `tol·up` floor at eps2 (ModularFunctionMinimizer) | 530–534 | ✓ |
| n==0 / seed-invalid / edm<0 gates (77–92) | 547–582 | ✓ (relaxed seed gate) |
| do-while outer; call inner (111–118) | `while iterate` 690 | ✓ inlined |
| Strategy≥1 HESSE `S==2 ‖ (S==1 && Dcovar>0.05)` (138–142) | 888–900 | ✓ |
| invalid Hessian → break (150–153) | 904–911 | ✓ |
| re-iterate if `edm>edmval && edm≥\|eps2·fval\|` (160–168) | 927–932 | ✓ exact |
| `maxfcn_eff = int(maxfcn·1.3)` on pass 0 (182–183) | 937–939 | ✓ |
| final `edm>10·edmval → MnAboveMaxEdm` (189–198) | 950, 952 | ✓ |
| call-limit `nfcn≥maxfcn → MnReachedCallLimit` (350–354) | 949 | ✓ |
| inner edm classification `<machine`/`<10·edmval`/else (356–368) | folded into `above_max` 950 | ✓ |

### Findings

- **Deliberate documented divergences (not bugs):**
  1. *Status-gated entry shortcut* (migrad.jl:720–722): skips the inner-loop
     body when `edm ≤ edmval && status == MnHesseValid`; C++ is a strict
     `do{...}while`. The load-bearing PR #10 / DAVIDON-audit subtlety — the
     shortcut fires *only* for an already-converged trustworthy-V warm restart
     (the MINOS/contour no-op case); for a placeholder-V seed (status ≠
     MnHesseValid) it does not fire, preserving do-while semantics (the IAM
     x_jm → 322 walk). Correctness-preserving optimization.
  2. *Relaxed seed-validity gate* (573–577): structural validity (params /
     gradient set, error available) vs C++'s effectively-no-op `seed.IsValid()`.
     More correct — accepts a bailed-but-usable `_hesse_diagonal_failure` seed.

- **MINOR (efficiency, not correctness): missing the C++ "2nd-pass invalid →
  bail" guard** (C++ 127–132: `if (ipass>0 && !min.IsValid()) return`). JuMinuit
  re-iterates under the same edm condition but lacks this early-out, so a
  non-converging fit **at Strategy ≥ 1** can run extra HESSE+DFP passes (bounded
  by the 1.3× call limit) before giving up, where C++ stops at pass 2. Same
  final verdict (invalid); JuMinuit spends more FCN calls. Narrow (S≥1
  non-converging only); ~3-line guard would restore exact parity.

- **Negligible:** at the no-improvement exit JuMinuit keeps `s0`'s old fval;
  C++ (size>1) records `pp.Y()` — differ by ≤ `eps·|fval|` (that branch's own
  entry condition), machine-precision.

- **Structural equivalences:** two-method split → one inlined loop; C++ `result`
  vector + reduced-state storage → JuMinuit `history` (storage-level-gated) +
  `final=s0`; MnPosDef bail returns a `FunctionMinimum` (C++) vs breaks-then-
  builds (JuMinuit).

- **Collaborators** (verified separately): `DavidonErrorUpdator`→davidon.jl and
  `VariableMetricEDMEstimator`→edm.jl line-by-line in `DAVIDON_CXX_AUDIT.md`;
  `MnLineSearch`+`MnParabola*`→linesearch.jl, `MnPosDef`→posdef.jl ported.

**Verdict: faithful port.** Every branch and exit path of both methods maps
correctly. Substantiates the `IAM_CONVERGENCE_GAP.md` § Fidelity claim
("core MIGRAD is faithful") with line-by-line evidence, consistent with the
Rosenbrock/Quad exact-match. Only non-cosmetic items: the deliberate
status-gated shortcut and the minor missing 2nd-pass-invalid bail.

---

## 3. MnMinos

`MnMinos.cxx` (213 lines) sets up each ±σ scan and delegates the actual
root-finding to `MnFunctionCross.cxx` (512 lines). JuMinuit splits these the
same way: `src/minos.jl` (the `FindCrossValue` setup + MinosError assembly) and
`src/function_cross.jl::_cross_core` (the parabolic root-find, shared with
MnContours). `function_cross.jl` is larger (1597 lines) because it also serves
contours, multi-fixed-parameter scans, the AD path, and warm-restart reuse.

### 3a. MnMinos::FindCrossValue (C++ MnMinos.cxx:94–197 ↔ minos.jl `minos(...)`)

| C++ branch | JuMinuit | Verdict |
|---|---|---|
| `err = dir·Error(par)`, `val = value + err` (119–120) | `sigma_i = √(2·up·V[ii])` (226), dir applied in `function_cross` | ✓ |
| limit clamp of `val` (122–129) | bounded-path int↔ext clamp (275–302) | ✓ (+ hardening below) |
| `xunit = √(up/m(ind,ind))`; other-param pre-shift `xt(i)+dir·xunit·m(ind,i)` (140–165) | `shift = σ·V[ik]/V[ii]`, seed_upper/lower (271) | ✓ **algebraically verified** (the 2·up & 2× factors cancel; minos.jl:234–238) |
| `upar.Fix(par); SetValue(par,val)` (167–168) | par_idx is the fixed scan param in `function_cross` | ✓ |
| `MnFunctionCross(...)` (172–173) | `function_cross(fmin, cf, par_idx, ±1; …)` (333, 367) | ✓ |
| AtMaxFcn / NewMinimum / AtLimit / !IsValid warnings (178–192) | MnCross flags + invalid-side ±σ placeholder (341–350) | ✓ (matches `MinosError::Upper/Lower`) |
| `maxcalls==0 → 2·(nvar+1)·(200+100n+5n²)` (111–114) | high-level default `maxcalls=1000` (minuit.jl:846) | **✗ divergence (below)** |

### 3b. MnFunctionCross (C++ MnFunctionCross.cxx ↔ function_cross.jl `_cross_core` + helpers)

| C++ branch | JuMinuit | Verdict |
|---|---|---|
| `aim = aminsv+up`, `tlf = tlr·up`, `tla = tlr`, `maxitr=15` (45–50) | 242, 261, `tla_base`, `maxitr` | ✓ |
| inner `MnMigrad(…, MnStrategy(max(0,strategy−1)))` (106) | `Strategy(max(0, level−1))` (799, 965) | ✓ exact |
| 1st MIGRAD; `flsb[0]=max(Fval,aminsv+0.1·up)`; `aopt=√(up/(f−fmin))−1` (119–142) | 270–276 | ✓ |
| converged `\|flsb[0]−aim\|<tlf` (143–144); clamp `[−0.5,1]` (146–149) | 278–281 | ✓ |
| 2nd MIGRAD; `dfda=(f1−f0)/(a1−a0)` (164–184) | 284–302 | ✓ |
| L300 `dfda<0` extend `aopt=alsb[0]+0.2·(it+1)` (188–242) | `while dfda<0`, `a[1]+0.2·count` (312–335) | ✓ |
| L460 linear extrap `aopt=alsb[1]+(aim−flsb[1])/dfda`; converge `adist<tla && fdist<tlf`; `[bmin,bmax]` clamp (244–266) | 343–355 | ✓ |
| 3rd MIGRAD + 3-point `noless` dispatch (288–351) | 357–404 | ✓ (incl. the "new straight line" L460-reentry, review BLOCKING #2) |
| L500 parabola loop: `MnParabolaFactory` fit, solve `=aim`, positive-slope root, converge at `ibest`, window/bad-point mgmt, replace worst (353–503) | `_parabola_fit3`/`_parabola_solve_for_aim`/`_three_point_classify` + L500 `while ipt<maxitr` (406–503) | ✓ line-cited |
| exits CrossNewMin / CrossFcnLimit / CrossParLimit / invalid / converged | `new_min` / `fcn_limit` / `par_limit` / `valid=false` / `valid=true` | ✓ (par_limit structural, below) |

### Findings

- **✗ Divergence (MODERATE, drop-in-compat): default MINOS call budget.** C++
  (and iminuit) default `maxcalls=0` → `2·(nvar+1)·(200+100·nvar+5·nvar²)`
  (≈30 100 for n=9); JuMinuit's high-level `minos!`/`minos` default is a fixed
  `maxcalls=1000` (minuit.jl:846, minos.jl:200). On larger fits JuMinuit MINOS
  can hit `fcn_limit` where C++/iminuit would keep going. User-overridable via
  `maxcall=`. **Recommended fix**: when `maxcall==0`, compute the C++ n-scaled
  default instead of falling back to 1000 (~3 lines; restores drop-in parity).

- **Structural-but-equivalent: `par_limit`/`aulim` detection.** C++ computes
  `aulim` inside MnFunctionCross with inline per-probe `limset && Fval<aim →
  CrossParLimit` exits (66–104, 135, 178, 227, 294, 495). JuMinuit's core
  `_cross_core` is limit-agnostic (operates in the caller's frame); the bounded
  wrapper detects `par_limit` via the int↔ext transform + a post-hoc aulim-style
  check (function_cross.jl:1291, 1370–1388). Same outcome (par_limit raised when
  the crossing lies beyond a bound); the *timing* of detection within the loop
  differs. Documented (function_cross.jl:1165–1168).

- **Hardening beyond C++ (not a gap):** the other-parameter pre-shift adds a
  sin-transform saturation pre-clamp for doubly-bounded params (minos.jl:254–302)
  to prevent `sin()` aliasing on large pre-shifts — a safety branch C++ lacks.

- **Extension beyond C++ (not a gap):** `sigma=k` k-σ MINOS errors (the
  `aopt·σ_i` scaling); C++ `MnMinos` is 1σ-only.

**Verdict: faithful port.** The root-finding core (`_cross_core`) is a
meticulous, C++-line-cited reproduction of MnFunctionCross — every branch
(L300/L460/L500, the noless dispatch, parabola fit, window/bad-point management)
and every exit (new-min / call-limit / par-limit / invalid / converged) maps,
with the inner-MIGRAD `Strategy−1` reduction and the covariance cross-correlation
pre-shift algebraically verified. The only substantive divergence is the
**smaller default call budget** (1000 vs n-scaled) — a drop-in-compat concern on
larger fits, easily fixed.

---

## Summary across the three audits

| Algorithm | Verdict | Substantive items |
|---|---|---|
| **MnHesse** | faithful | bounded-parameter step clamp not implemented (`has_limits=false`; documented Phase-1 deferral; unbounded fits unaffected) |
| **MIGRAD** | faithful | deliberate status-gated entry shortcut (correctness-preserving); missing C++ 2nd-pass-invalid early-bail (efficiency-only, S≥1 non-converging) |
| **MnMinos** | faithful | default call budget 1000 vs C++/iminuit n-scaled (drop-in-compat; ~3-line fix) |

All three are faithful ports of the C++ Minuit2 algorithm. No whole branch is
silently absent; the divergences are (a) documented deliberate optimizations,
(b) a documented Phase-1 bounds deferral, (c) a narrow efficiency gap, and (d) a
default-budget mismatch. Two carry a concrete, small recommended fix (MnHesse
bounded clamp; MnMinos default budget); one is a deliberate keep (MIGRAD
shortcut); the rest are same-result reformulations or hardening beyond C++.
