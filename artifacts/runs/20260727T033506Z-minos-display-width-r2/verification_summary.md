# MINOS fixed-width display verification record

- Run ID: `20260727T033506Z-minos-display-width-r2`
- Branch: `codex/minos-limit-display-fix`
- Baseline commit: `1c845b27f928aea8fe66d91a41c0ae794336c884`
- Status: verified

## Change

The single-result `text/plain` MINOS box now:

- labels an invalid side as `invalid` in the shared summary row, avoiding
  redundant `upper invalid` and `lower invalid` prefixes;
- preserves the explicit `Upper/Lower at limit` and `distance` wording in
  the detail row;
- uses a compact scientific representation only when an at-limit distance
  would exceed the nine-character numeric budget of a 35-column cell.

The stored `MinosError` values, validity flags, closed-interval semantics,
compact display, HTML, LaTeX, plotting, serialization, and MINOS numerical
algorithms are unchanged.

## Verification

| Check | Result | Evidence |
| --- | --- | --- |
| Pre-fix reproduction | 118 passed, 3 failed; the three new width cases failed | `pre_fix_failing.log` |
| Focused MINOS presentation tests | 121/121 passed | `focused.log` |
| Rich display tests | 151/151 passed | `display.log` |
| Rendered-box inspection | All nine nonempty lines in each of three representative boxes have text width 73 | `presentation_check.log` |
| At-limit width sweep | 30,033 upper-, lower-, and double-limit boxes passed across special values and randomized magnitudes | `width_sweep.log` |
| Documentation build | Completed successfully; existing warning-only classes remain | `docs_build.log` |
| Full package tests | 5030 passed, 1 expected broken, 0 failed, 0 errored | `pkg_test.log` |
| Diff hygiene | `git diff --check` passed | command output at closeout |

The full suite includes 99/99 core MINOS tests, 62/62 C++ JSON MINOS oracle
tests, 19/19 bounded Gaussian C++ MINOS oracle tests, and 121/121 public
at-limit presentation tests.

## Pre-commit research-integrity check

- M1: reproduced three concrete failures before the fix; tested mixed
  validity, both-invalid, very large distances, and a 30,033-case
  magnitude sweep. The full test suite and C++ oracle suites reject the
  counter-hypothesis that numerical MINOS behavior changed.
- M2: not applicable; no citations were added or changed.
- M3: not applicable; no measured or literature-derived numerical result
  was introduced. The 73-column layout invariant was measured directly.
- M4: not applicable; no literature relationship claim was made.
- M5: checked by inspecting the diff and running the MINOS oracle and full
  package suites. The change is confined to presentation strings and a
  presentation-only formatter; no scientific observable or plot data
  changed.
- M6: every verification claim above is bound to an executed command and
  its persisted log in this directory.
- M7: tested the opposing framing that this was a general arbitrary-field
  box-truncation problem. The reproduced defect is limited to the new
  at-limit/invalid status strings, so the fix and its claims remain scoped
  to those rows rather than changing shared box formatting globally.
- M8: not applicable; no production computation or reusable numerical
  method was started.

No Nullius approval was pending at session recovery, so no
`integrity-record` approval receipt was required.

## SHA-256

```text
5c177b9d89fb671108d9b693b4099ef66417338988c0805bee326cdf941cfe55  display.log
f51102c516d7b172c6460949c76f307da7eb7d7b089cf8674d530bb156deabf0  docs_build.log
63c830c894b6d42c6bbf4f08851cc3f34ebfafaa9955da7091750d271e9e5b60  focused.log
23e02335da4b810cd8c78b63a04d8dd7b3b81beff3c2d6d7da74eccb4de404cc  pkg_test.log
9159fd9af7eef553455d3b63caea4f2e4887f6c5605c5af3e1c47a7ef4337b20  pre_fix_failing.log
49fa502706c15b592474e7b000a25517f0da0620a175fca483bc69fcfaf4abe7  presentation_check.log
6c288cc96d1774d12ad48d87f9e853f478189def41a90239ee1bf47636876ea8  width_sweep.log
```
