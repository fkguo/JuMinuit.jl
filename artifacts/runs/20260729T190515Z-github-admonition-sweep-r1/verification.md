# Portable admonition rendering verification

The repository-wide Markdown audit found 13 active Documenter-style note or
warning blocks across eight public documentation pages. Each was converted to a
standard Markdown blockquote with an explicit bold `Note —` or `Warning —`
title. The wording and scientific content were not changed.

Checks:

- A repository-wide anchored scan found no remaining active line-start
  `!!! note` or `!!! warning` syntax in Markdown files.
- GitHub's GFM rendering API was applied independently to all eight changed
  pages. It rendered all 14 portable note/warning titles (13 converted here and
  one converted previously) inside blockquotes and rendered no literal
  Documenter admonition syntax.
- `julia --project=docs docs/make.jl` completed successfully.
- The generated Documenter HTML contains all 14 portable note/warning titles
  and no literal Documenter admonition syntax.
- `git diff --check` passed.

The Documenter build retained pre-existing, unrelated warnings about unresolved
API references, docstrings not included in the manual, generated page size, and
SVG fallbacks.
