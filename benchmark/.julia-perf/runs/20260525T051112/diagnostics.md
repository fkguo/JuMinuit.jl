# Julia Perf Diagnostics

- timestamp: 2026-05-25T05:11:12Z
- mode: standalone
- verdict: FAIL

## Hard Failures
- missing-baseline: No baseline file and --save-baseline not provided

## Soft Warnings
- suite-warning: failed to include suite file: LoadError: ArgumentError: Package JuMinuit not found in current path.
- Run `import Pkg; Pkg.add("JuMinuit")` to install the JuMinuit package.
in expression starting at /Users/fkg/Coding/Agents/ResearchWork/JuMinuit/benchmark/bench_migrad_suite.jl:26

## Notes
- No benchmark rows produced
