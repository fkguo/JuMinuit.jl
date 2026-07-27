# Research-integrity closeout

## M1 — implementation bug passing self-review

Checked the final diff by function and searched every call site of
`_format_value_minos`, `_format_minos_err`, `_minos_show_side`, and
`_minos_plot_error`.

Counter-hypothesis: a boundary displacement might still enter the ordinary
asymmetric formatter through the internal compact `_value_cell` helper even
though the main fit table uses `_split_cells`. This path was possible for a
direct helper call, so the helper was hardened to reject either at-limit side,
and a regression test now verifies that it falls back to the Hesse display.
The final related and full suites pass.

The test-first negative control also rejected the old implementation before
the display changes, demonstrating that the new tests do not pass vacuously.

## M2 — hallucinated citation

Not applicable: no citation or bibliographic claim was added.

## M3 — hallucinated measurement or result

Not applicable to external measurements. All reported test counts are taken
from the persisted command logs. The real-fit regression values are computed
inside the tests and are not promoted as scientific measurements.

## M4 — shortcut reliance

Not applicable: no literature-relationship claim was made.

## M5 — bug as insight

The opposing hypothesis was that the observed at-limit number might be a
closed confidence error. Real bounded profiles disconfirm this: wide bounds
produce two crossings and a closed interval, while flat or one-sided bounded
profiles retain valid computation flags but set the independent parameter-limit
flag and do not close the interval. Full MINOS and C++ oracle tests pass.

The algorithm and serialization diffs are empty, so no changed numeric result
is being interpreted as a new finding.

## M6 — methodology fabrication

The exact related-test, full-suite, and documentation-build commands and their
outputs are persisted in this run directory. No unexecuted verification is
claimed.

## M7 — frame lock

Reframed the issue from the opposite interpretation: “the best-fit parameter
is at the boundary.” The flat-profile regression has a strictly interior
best-fit value while both profile scans reach their limits, so that framing is
false. The final displays therefore report best-fit proximity and MINOS scan
termination as separate events.

## M8 — reinvention over available reuse

Not applicable: this change adds no numerical method, solver, or production
computation. It reuses the existing MINOS status fields and existing test
suite.
