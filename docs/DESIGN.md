# Design notes — JuMinuit.jl

Architectural decisions, written down as they're made so the *why* is
preserved. Each entry is small (1–10 lines).

This is a living document. Phase 0 will populate the early sections;
later phases append.

> Format: `## DR-NNN: short title — DECIDED YYYY-MM-DD`
> Then: **Context** | **Decision** | **Consequences** | **Revisit when**.

---

## DR-001: Repository layout — DECIDED 2026-05-25

**Context**: standard Julia-package layout vs. monorepo with multiple
sub-packages.

**Decision**: single package, standard layout — `src/`, `test/`,
`docs/`, `benchmark/`, `tools/`, `scripts/`. C++ reference at
`reference/Minuit2_cpp/` (gitignored, mirror only). Skill snapshot at
`.claude/skills/julia-perf/` (vendored from autoresearch-lab).

**Consequences**: simple to maintain; one registry entry when Phase 0
ships. If we later want a separate `JuMinuitCompat.jl` for IMinuit-
compatibility shim, it ships as a submodule, not a separate package
(at least until Phase 3 if user demand surfaces).

**Revisit when**: Phase 3 IMinuit-compat work is in flight and the
submodule approach feels too tight.

---

## DR-002: Pin to upstream Minuit2 v6.24.0 — DECIDED 2026-05-25

**Context**: how to handle upstream C++ Minuit2 evolution. New versions
may shift bits at the 1e-13 ULP and break reference data.

**Decision**: pin to GooFit/Minuit2 commit `57dc936` (v6.24.0). Upstream
bumps happen at major Julia version boundaries with a checklist in
`tools/regen_reference.md`.

**Consequences**: reproducibility wins; we lock in any C++ bugs
existing at 6.24.0; cross-platform IEEE drift addressed via the
tolerance hierarchy in `tools/regen_reference.md`.

**Revisit when**: JuMinuit v2.0 (planned upstream bump), or if a
critical-correctness bug surfaces upstream that we want to inherit.

---

## DR-003: Strict mirror of C++ bounds transforms (sin/sqrt) — DECIDED 2026-05-25

**Context**: ROADMAP § 10 Q2. Modern alternatives (tanh, logistic)
exist; sin/sqrt have known pathology near bounds.

**Decision**: strict mirror for v1.0. Alternatives in Phase 2 if
demanded.

**Consequences**: iminuit copy-paste compatibility preserved; users
hitting the sin-pathology can use unbounded parameters with a manual
penalty.

**Revisit when**: a user submits a minimal repro of the sin-pathology
biting a real HEP fit *and* demonstrates tanh would have fixed it.

---

## DR-004: Dense `Symmetric{Float64,Matrix{Float64}}` as default — DECIDED 2026-05-25

**Context**: ROADMAP § 10 Q3. C++ uses packed lower-triangular storage
for the inverse Hessian; Julia can do either.

**Decision**: dense `Symmetric{Float64,Matrix{Float64}}` is the Phase 0
default. Packed variant available behind a feature flag in `linalg.jl`
for benchmark comparison.

**Consequences**: OpenBLAS DSYMV vectorizes better than DSPMV at the
typical n ≤ 50; both representations fit in L1. The factor-of-2 memory
cost is irrelevant. `BLAS.syr!` requires `parent(::Symmetric)` to get
the underlying `Matrix` — document the upper/lower triangle convention
in `linalg.jl`.

**Revisit when**: Day-26 benchmark shows packed beating dense by > 10%
on any blocking scenario.

---

## DR-005: Threaded gradient opt-in, not default — DECIDED 2026-05-25

**Context**: ROADMAP § 10 Q5. iminuit and C++ default to single-threaded
(OpenMP opt-in via CMake flag).

**Decision**: `migrad(...; threaded_grad=false)` default. Phase 2.2 adds
`threaded_grad=true`; requires a *pure* FCN.

**Consequences**: predictable for benchmarking parity; matches iminuit
default; users who know they have pure expensive FCNs opt in. BLAS
thread count managed independently.

---

## DR-006: `@enum CovStatus` for error-matrix status — DECIDED 2026-05-25

**Context**: ROADMAP § 10 Q7. C++ uses tag types; Julia options were
`@enum`, `Symbol`, `Val{:tag}`.

**Decision**: `@enum CovStatus { MnHesseFailed, MnMadePosDef,
MnInvertFailed, MnNotPosDef, MnHesseValid }`.

**Consequences**: inferable, printable, dispatchable on. Slight
overhead vs `Symbol` is negligible; the inferability win in switch-on-
status code paths matters.

---

## DR-007: API parity — both `m.values` and `values(m)` — DECIDED 2026-05-25

**Context**: ROADMAP § 10 Q6. iminuit users expect `m.values`; Julia
idiom is `values(m)`.

**Decision**: ship both. `Base.getproperty(m, :values)` returns
`values(m)` internally; the function-style is the canonical form.

**Consequences**: iminuit/IMinuit.jl copy-paste works; idiomatic Julia
also works; one-line cost per accessor.

---

## DR-008: Phase 0 locks Strategy = 0 — DECIDED 2026-05-25

**Context**: C++ `VariableMetricBuilder.cxx:138–173` invokes `MnHesse`
internally when `Strategy ≥ 1 && Dcovar > 0.05`. iminuit default is
Strategy 1. Phase 0 doesn't have HESSE yet.

**Decision**: Phase 0 ships Strategy = 0 only. Users who request
Strategy ≥ 1 get an `ArgumentError("Strategy ≥ 1 requires Phase 1
(HESSE). Use Strategy(0).")`.

**Consequences**: avoids a missing-code-path bug; documents the gate
between Phase 0 and Phase 1 cleanly. Strategy 0 corresponds to "fast
mode" in iminuit — a defensible Phase 0 minimum.

**Revisit when**: Phase 1 ships `hesse.jl` and the inner-HESSE call
path.

---

## DR-009: `storage_level=0` default in Phase 0 — DECIDED 2026-05-25

**Context**: C++ default is `storage_level=1` — appends a `MinimumState`
per iteration to `FunctionMinimum::history`. Per-iteration allocation
breaks the §3.4 zero-alloc gate if it's the default.

**Decision**: Phase 0 defaults `storage_level=0` (seed + final state
only). Tests requiring full trace set it explicitly. Phase 1 exposes
the keyword to users.

**Consequences**: zero-alloc gate stays achievable; users wanting full
trace explicit about it; matches the spirit of the C++ default
(history) without paying the cost on every run.

---

## (Open) DR-010: License — LGPL 2.1 vs MIT

**Context**: ROADMAP § 10 Q13. C++ Minuit2 is LGPL 2.1. A literal line-
by-line port creates derivative work; Julia ecosystem default is MIT.

**Status**: open. Requires user decision before src/ code lands.

**Options under consideration**:
- (a) LGPL 2.1 for JuMinuit.jl (preserves derivative-work rules)
- (b) Clean-room from FORTRAN Minuit paper (James 1975) + iminuit Python
  source only, NOT reading the C++ — but the C++ remains the only
  complete numerical oracle
- (c) Mixed: oracle data via C++ build (gitignored, not redistributed),
  Julia source written from algorithmic descriptions

**Blocks**: Phase 0 src/ code. Reference data generation is unblocked
since reference data is computed output, not source.

---

## Future decisions (placeholders — to be made in Phase 0/1)

- DR-011: BLAS vendor for the gate (OpenBLAS / MKL / Accelerate)
- DR-012: Designated reference machine
- DR-013: Minimum Julia version + `[compat]` policy
- DR-014: AD backend primary (ForwardDiff vs Enzyme) — Phase 2.1 lock
- DR-015: SVector/MVector small-n specialization — Phase 0 day-26
  decision
