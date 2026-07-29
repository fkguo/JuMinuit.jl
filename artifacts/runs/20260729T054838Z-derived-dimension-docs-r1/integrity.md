# Research-integrity check

- Run ID: `20260729T054838Z-derived-dimension-docs-r1`
- Boundary: durable NativeMinuit documentation and API docstring update
- Pending Nullius approval at recovery: none

## M1 — implementation or reasoning passing self-review

The load-bearing counter-hypothesis was that the number of displayed outputs
always equals the likelihood-ratio degrees of freedom. The two independent
derivations refuted that simplification: the relevant quantity is the local
rank on identifiable fit directions. The final text states the rank condition,
including locally dependent outputs, and does not change executable code.

A repository-wide search found no remaining statement that `ndof` is simply
the total number of fit parameters. The source diff in
`src/error_sampling.jl` is confined to the `delta_chisq` docstring.

## M2 — hallucinated citation

Not applicable. No new paper, author, identifier, quotation, or bibliographic
claim was added. The new guide paragraph links to the existing NativeMinuit
contour-convention tutorial and its existing references.

## M3 — hallucinated numerical result

The three displayed numerical statements were evaluated with the current
NativeMinuit implementation:

- `delta_chisq(0.6826894921370859, 1) = 1.0000000000000002`;
- `delta_chisq(0.6826894921370859, 2) = 2.295748928898635`;
- `chisq_cl(1.0, 2) = 0.3934693402873665`.

The executed output is in `numerical_thresholds.log`; the threshold and inverse
functions also passed all 193 tests in `test_error_sampling.jl`.

## M4 — unsupported literature relationship

Not applicable. No citation-network, precedence, influence, or paper-to-paper
relationship claim was introduced.

## M5 — numerical or implementation artifact presented as insight

No numerical algorithm or statistical implementation was changed. The
threshold values were checked through both the analytic derivations and the
package implementation. The complete targeted error-sampling test file passed
193/193, and the documentation built from the current local checkout.

## M6 — methodology fabrication

Every verification claim is bound to an executed artifact:

- independent derivations and comparison:
  `derivation_summary.md` and `derivation_verification_matrix.json`;
- direct threshold calculation: `numerical_thresholds.log`;
- targeted package tests: `test_error_sampling.log`;
- local-source Documenter build: `docs_build.log`;
- rendered-text and cross-link checks: `doc_semantics_check.log`.

The deterministic Markdown checker reported many pre-existing, intentional
code-span-versus-math style candidates in the two long documentation files.
No bulk fixer was applied. The actual target renderer, Documenter, completed
successfully and rendered the new mathematics and cross-link correctly.

## M7 — frame lock

The opposing framing was tested explicitly: neither the total number of fit
parameters nor the raw number of output values or curve-grid points determines
the threshold in general. The final wording uses the effective local rank,
distinguishes separate scalar intervals from a joint two-dimensional pole
contour, and distinguishes pointwise from simultaneous curve coverage.

## M8 — reuse before heavy computation

Not applicable. No heavy or production computation and no reusable numerical
method were started; this change is documentation-only.
