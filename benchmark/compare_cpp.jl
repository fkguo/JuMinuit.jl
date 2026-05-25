# SPDX-License-Identifier: LGPL-2.1-or-later
#
# Julia ↔ C++ Minuit2 wall-time comparison.
#
# Phase 0 §3.4 Criterion 2 driver. Reads:
# - Julia benchmark medians from
#   `benchmark/.julia-perf/runs/latest/benchmarks.json`
#   (produced by `scripts/run_gate.sh`).
# - C++ benchmark medians from the cpp_bench JSON output
#   (`benchmark/cpp/build/cpp_bench`).
#
# Emits a ratio table and a verdict:
# - PASS if every Julia / C++ ratio ≤ 1.5 (Phase 0 gate).
# - WARN if ratio is 1.5-1.6 (within instrumentation noise).
# - FAIL if any ratio > 1.6.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "scripts"))
using JSON3
using Printf

const JULIA_PERF_DIR = joinpath(@__DIR__, ".julia-perf", "runs", "latest")
const CPP_BENCH = joinpath(@__DIR__, "cpp", "build", "cpp_bench")

function load_julia_results()
    path = joinpath(JULIA_PERF_DIR, "benchmarks.json")
    isfile(path) || error("Julia benchmarks.json not found at $path. " *
                          "Run `scripts/run_gate.sh --save-baseline` first.")
    data = JSON3.read(read(path, String))
    return Dict(String(b.name) => Float64(b.current_median_ns) for b in data)
end

function load_cpp_results()
    isfile(CPP_BENCH) ||
        error("$CPP_BENCH not built. Run `cmake --build benchmark/cpp/build` first.")
    output = read(`$CPP_BENCH`, String)
    data = JSON3.read(output)
    return Dict(String(b.name) => Float64(b.median_ns) for b in data)
end

function main()
    julia_t = load_julia_results()
    cpp_t = load_cpp_results()
    shared = intersect(Set(keys(julia_t)), Set(keys(cpp_t)))

    println("\nPhase 0 §3.4 Criterion 2 — Julia vs C++ Minuit2 wall time")
    println("Median per migrad() call (BenchmarkTools median; n_samples=50)")
    println("Strategy(0); BLAS.set_num_threads(1) on both sides; Apple M3.")
    println()
    println(rpad("Benchmark", 22), rpad("Julia (μs)", 14),
            rpad("C++ (μs)", 14), rpad("Julia/C++", 12), "Gate (≤1.5×)")
    println("-"^76)

    max_ratio = 0.0
    fail_count = 0
    warn_count = 0
    for name in sort(collect(shared))
        j = julia_t[name] / 1000   # μs
        c = cpp_t[name] / 1000
        ratio = j / c
        max_ratio = max(max_ratio, ratio)
        status = if ratio <= 1.5
            "✓"
        elseif ratio <= 1.6
            warn_count += 1
            "⚠"
        else
            fail_count += 1
            "✗"
        end
        @printf "%-22s%-14.2f%-14.2f%-12.3f%s\n" name j c ratio status
    end
    println("-"^76)
    @printf "Maximum ratio: %.3f\n" max_ratio

    verdict = if fail_count > 0
        "FAIL"
    elseif warn_count > 0
        "WARN"
    else
        "PASS"
    end
    @printf "\nVerdict: %s\n" verdict
    return verdict == "PASS" ? 0 : (verdict == "WARN" ? 2 : 1)
end

exit(main())
