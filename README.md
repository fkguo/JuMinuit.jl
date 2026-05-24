# JuMinuit.jl

Native Julia port of the C++ [Minuit2](https://root.cern.ch/doc/master/Minuit2Page.html) function minimization library, targeting C++-comparable performance and a Julia-native API.

**Status**: Phase 0 (proof of concept). See [ROADMAP.md](ROADMAP.md).

## Goals

- Drop the C++ dependency that [IMinuit.jl](https://github.com/fkguo/IMinuit.jl) and [iminuit](https://github.com/scikit-hep/iminuit) wrap.
- Match GooFit/Minuit2 standalone performance on representative HEP workloads (≤ 1.5× C++ time, Phase 0 gate).
- Add Julia-native features unavailable in C++: AD-driven gradients (ForwardDiff/Enzyme), Threads-parallel likelihoods, recipes for plotting.

## Reference

The C++ source at [GooFit/Minuit2](https://github.com/GooFit/Minuit2) is mirrored locally at `reference/Minuit2_cpp/` (gitignored) for porting reference and benchmark comparisons.

## Layout

```
src/         Julia source (TBD per ROADMAP Phase 0)
test/        Julia test suite
docs/        Design notes, porting tables, benchmark history
reference/   C++ Minuit2 mirror (not tracked)
```
