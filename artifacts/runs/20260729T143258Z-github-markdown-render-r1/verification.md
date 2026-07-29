# GitHub and Documenter rendering verification

The Documenter-only `!!! note` syntax was parsed as ordinary text followed by
an indented code block on GitHub. The note now uses a standard Markdown
blockquote, which both renderers support.

The GitHub GFM endpoint recognized:

- the note block;
- all three list items;
- all nine inline mathematical expressions as `math-renderer` elements;
- the relative documentation link without interpreting `@ref` as a username.

No raw `$...$` expression remained outside a math-renderer element in the
target block.

The local Documenter build completed successfully. Its generated HTML contains
the same note block, three list items, nine mathematical spans, and a rewritten
link to `tutorials/minos_contours.html`. The build emitted only the existing
API-reference, docstring-coverage, page-size, and SVG-fallback warnings.

The standalone Markdown hygiene checker still reports 113 pre-existing
whole-file findings, mostly intentional code spans and Documenter-specific
cross-references elsewhere in the page. None is introduced by the corrected
note. The scientific wording and numerical statements are unchanged.
