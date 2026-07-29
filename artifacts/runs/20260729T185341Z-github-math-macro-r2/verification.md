# GitHub math-macro verification

GitHub recognized the original `$...$` spans but rejected `\operatorname` during
browser-side typesetting. The formulas now use `\mathrm{Re}` and
`\mathrm{Im}`, which preserve the intended notation in both GitHub MathJax and
Documenter KaTeX.

The exact GitHub page at commit `d3f86e8` was inspected after the push. The
derived-target note contains three list items and nine mathematical
expressions. All nine expressions produced MathML nodes. Neither the note nor
the rest of the page reported a disallowed-macro error.

The GFM parser independently found no unrendered `$...$` expression in the
note, the later complex-pole support-function paragraph, or the corresponding
NativeMinuit usage-skill paragraph. A source scan found no remaining
`\operatorname` in those two Markdown files.

The local Documenter build also completed successfully with only the existing
API-reference, docstring-coverage, page-size, and SVG-fallback warnings. The
scientific content and numerical thresholds were unchanged.
