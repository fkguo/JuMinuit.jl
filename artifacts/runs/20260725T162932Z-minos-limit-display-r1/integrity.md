# MINOS at-limit presentation verification

## Scope and invariant

This change is limited to public presentation and export semantics. It does
not alter the MINOS search, bound-distance calculation, termination flags,
state snapshots, or C++/iminuit numerical compatibility.

For an at-limit side, the raw `upper` or `lower` field remains the signed
distance from the minimum to the parameter bound. Public displays label that
number as a distance to the limit, and plot recipes omit it from confidence
whiskers. A closed interval requires two finite $\Delta\mathrm{FCN}$
crossings.

## Test-first evidence

The initial test-first execution against the unmodified presentation code
reported 58 passes, 30 failures, and 20 errors. The failures/errors localized
the missing behavior to the public helper, text/HTML/LaTeX presentation,
plot whiskers, and serialized parameter-limit flags. This interactive
baseline output was not retained as a separate log.

The final focused test covers upper-only, lower-only, both-bounds, a normal
internal crossing, and one invalid side. It also checks the exact motivating
regression value: a minimum of `3.7799991402523414` below an upper bound of
`3.78`.

## Verification results

- Focused presentation test: 112/112 passed; see
  [`targeted_final.log`](targeted_final.log).
- Full `Pkg.test()`: 5004 passed, one repository-declared broken test, no
  failure or error; see [`pkg_test_offline.log`](pkg_test_offline.log).
- The full suite includes the existing unbounded C++ MIGRAD oracle (166/166),
  bounded C++ MIGRAD oracle (46/46), C++ MINOS JSON oracle (62/62), and bounded
  Gaussian MINOS oracle (19/19).
- Worktree Documenter build completed successfully; the new
  `has_closed_interval` API and tutorial references resolved. Remaining
  warnings are pre-existing warn-only references and page-size notices; see
  [`docs_build_final.log`](docs_build_final.log).
- `git diff --check` passed.

The first two online `Pkg.test()` attempts stopped before test execution
because both the Julia package server and GitHub clone endpoint ended their
TLS handshakes early. All test-target dependencies were already installed,
so the complete standard `Pkg.test()` was rerun with `Pkg.offline(true)`;
Julia resolved the full test target from the local cache and executed it.

## Research-integrity M1–M7

- **M1 — Counter-hypothesis search:** all `MinosError` consumers under `src/`
  and `ext/` were enumerated. The candidate counter-hypothesis was that a
  remaining public consumer could still treat a bound distance as a
  confidence error. The public display, export, plot, and serialization
  consumers are now status-aware. Numerical contour and MINOS code continues
  to consume the unchanged raw values.
- **M2 — Source support:** no new literature or scientific numerical claim is
  introduced. The preserved semantics are checked by the repository's
  existing C++ oracle tests and bounded MINOS tests.
- **M3 — Representativeness:** the focused cases span each side separately,
  both sides at limits, a genuine two-sided crossing, and a one-sided invalid
  result. Fit-table text/HTML, single-result text, LaTeX, single/vector plot
  recipes, and serialization are exercised.
- **M4 — Synthesis fidelity:** documentation distinguishes clean termination
  from a closed interval and does not reinterpret an at-limit distance as a
  small uncertainty.
- **M5 — Falsification and negative controls:** normal internal crossings
  retain their prior asymmetric format and whiskers. At-limit and invalid
  sides produce zero plot-whisker length, while their raw stored values remain
  unchanged. Legacy serialized payloads without the new keys remain readable.
- **M6 — Execution provenance:** reported test counts and warning status come
  from the retained logs in this run directory. The full suite and focused
  suite were executed in this worktree on branch
  `codex/minos-limit-display-fix`.
- **M7 — Opposing interpretation:** `is_valid(e) == true` was explicitly tested
  for at-limit results while `has_closed_interval(e) == false`. This prevents
  the plausible but incorrect inference that clean termination guarantees two
  $\Delta\mathrm{FCN}$ crossings.

No A1–A5 approval is requested or implied by this record.
