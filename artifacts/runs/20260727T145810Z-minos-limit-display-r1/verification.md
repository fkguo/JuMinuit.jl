# MINOS limit-display verification

## Scope

- Runtime display changes are limited to `src/display.jl`, the pretty-print
  section of `src/minos.jl`, the rich-output section of `src/minuit.jl`, and
  `src/plot_recipes.jl`.
- `src/serialize.jl`, `src/function_cross.jl`, `src/migrad.jl`, and
  `src/migrad_bounded.jl` have no diff from `main`.
- The existing `MinosError` fields, `is_valid`, `has_closed_interval`, MINOS
  computation, and returned numeric values are unchanged.

## Regression coverage

`test/test_minos_limit_display.jl` exercises real fits for:

1. A wide-bound quadratic with both requested confidence crossings found.
2. An interior best fit with a flat profile whose scans reach both limits.
3. A one-sided truncation with the opposite crossing preserved as a MINOS
   error.
4. A best-fit value close to a limit, with independent best-fit proximity and
   MINOS scan-termination warnings.

It also checks synthetic invalid and at-limit combinations, serialization
round-trips, fixed-width text output, HTML and LaTeX semantics, and plot
whisker omission.

## Test results

- Test-first negative control: the new regressions failed against the old
  wording and are recorded in `test-before-implementation.log`.
- Related final tests:
  `julia --project=. -e 'using LinearAlgebra, NativeMinuit, RecipesBase, Test; ...'`
  passed 627 tests. See `related-tests-post-integrity.log`.
- Final full suite:
  `julia --project=. -e 'using Pkg; Pkg.test()'`
  passed 5092 tests with 1 existing broken test. See `pkg-test-final.log`.
- Final documentation build:
  `julia --project=docs docs/make.jl`
  completed successfully. See `docs-build-final.log`. The log retains the
  repository's existing unresolved-reference, uncatalogued-docstring, and
  generated-page-size warnings.

## Static checks

- `git diff --check`: clean.
- No runtime display contains the old ambiguous phrases `upper at limit`,
  `lower at limit`, or `No params at limit`.
- Plot recipes return zero whisker length for every at-limit side and label
  the interval as truncated by a parameter limit.
