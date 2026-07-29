# GitHub and Documenter rendering verification

## Correction

This run verified GitHub's GFM parsing endpoint and the local Documenter output,
but it did not verify GitHub's subsequent browser-side macro whitelist. A live
GitHub rendering showed that `\operatorname` was rejected even though the GFM
endpoint had wrapped the expression in a `math-renderer` element. This run is
therefore superseded by
`20260729T185341Z-github-math-macro-r2`.

The Documenter-only `!!! note` syntax was parsed as ordinary text followed by
an indented code block on GitHub. The note now uses a standard Markdown
blockquote, which both renderers support.

The GitHub GFM endpoint recognized:

- the note block;
- all three list items;
- all nine inline mathematical expressions as `math-renderer` elements;
- the relative documentation link without interpreting `@ref` as a username.

No raw `$...$` expression remained outside a math-renderer element in the
target block at the parsing stage. That condition was necessary but not
sufficient for browser-side typesetting.

The local Documenter build completed successfully. Its generated HTML contains
the same note block, three list items, nine mathematical spans, and a rewritten
link to `tutorials/minos_contours.html`. The build emitted only the existing
API-reference, docstring-coverage, page-size, and SVG-fallback warnings.

The standalone Markdown hygiene checker still reports 113 pre-existing
whole-file findings, mostly intentional code spans and Documenter-specific
cross-references elsewhere in the page. None is introduced by the corrected
note. The scientific wording and numerical statements are unchanged.
