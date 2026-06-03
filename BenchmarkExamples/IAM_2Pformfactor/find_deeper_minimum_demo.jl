# SPDX-License-Identifier: LGPL-2.1-or-later
#
# `find_deeper_minimum` on the IAM ππ fit — WORKED EXAMPLE (demonstration, not a
# unit test).
#
# A single Strategy-1 MIGRAD from the published LEC start lands in a SHALLOW basin
# (χ² ≈ 379). `find_deeper_minimum`'s data-resampling dispatch climbs out of it:
# each round bootstrap-resamples the 85 phase-shift points, re-fits each resample
# (those drift toward whichever basin best explains that subset), clusters the
# results with `find_solution_modes(...; refine=true)`, adopts the deepest valid
# basin found on the ORIGINAL data, and repeats — here 379 → 285 → 279 → 260 → 255
# over four adopt-rounds (Δχ² ≈ 124, seed=1). Same mechanism as PHASE 1 of
# `error_crosscheck.jl` (which reaches ≈235 from a multi-start seed); from this
# colder single-start it lands in a different, comparably deep basin — exactly the
# path-sensitivity the heuristic warns about. ONE call replaces the hand loop.
#
# Companion: `error_crosscheck.jl` (the full multi-basin error-analysis study).
#
# Run (needs CSV/DataFrames/StaticArrays/QuadGK + JuMinuit):
#     julia --project=. BenchmarkExamples/IAM_2Pformfactor/find_deeper_minimum_demo.jl
#   Tunable env (defaults): IAM_NDISC=20  IAM_SEED=1
#
# Physics: GKPY Roy-equation ππ phase shifts; IAM with SU(3) NLO LECs.

using LinearAlgebra, Random, Statistics
BLAS.set_num_threads(1)
const NDISC = parse(Int, get(ENV, "IAM_NDISC", "20"))   # bootstrap resamples / round
const SEED  = parse(Int, get(ENV, "IAM_SEED", "1"))

const IAM_DIR = @__DIR__
cd(IAM_DIR)
using CSV, DataFrames, StaticArrays, QuadGK, JuMinuit

# ── constants + model (mirror error_crosscheck.jl / bench.jl setup) ──────────
const unit = 1.0
const fpi = 92.21unit; const mpic = 139.57018unit; const mpi0 = 134.9766unit
const meta = 547.862unit; const mkc = 493.677unit; const mk0 = 497.614unit
const mpi = (2mpic + mpi0)/3; const mk = (mkc + mk0)/2; const μ = 770.0unit; const ϵ = eps()
struct TwoBodyChannel{T<:AbstractFloat}; m1::T; m2::T; end
qon(s, m1, m2) = sqrt((s - (m1+m2)^2) * (s - (m1-m2)^2))/(2sqrt(s))
const ππ = TwoBodyChannel(mpi, mpi); const KK = TwoBodyChannel(mk, mk)
const ηη = TwoBodyChannel(meta, meta); const πη = TwoBodyChannel(mpi, meta)
const Kπ = TwoBodyChannel(mk, mpi); const Kη = TwoBodyChannel(mk, meta)
include(joinpath(IAM_DIR, "src", "init_const.jl"))
include(joinpath(IAM_DIR, "src", "amplitudes.jl"))
include(joinpath(IAM_DIR, "src", "tmatrix.jl"))
include(joinpath(IAM_DIR, "src", "unitarity_modification.jl"))
include(joinpath(IAM_DIR, "src", "phaseshifts.jl"))
const lecr0 = [0.56e-3, 1.21e-3, -2.79e-3, -0.36e-3, 1.4e-3, 0.07e-3, -0.44e-3, 0.78e-3]

_load(f) = JuMinuit.Data(DataFrame(CSV.File(f, header=[:w,:δ,:err], delim=' ', ignorerepeated=true)))
d00 = _load("./datajl/pipi/pipi00_Roy-GKPY_PRD83_074004.dat")
d11 = _load("./datajl/pipi/pipi11_Roy-GKPY_PRD83_074004.dat")
d20 = _load("./datajl/pipi/pipi20_Roy-GKPY_PRD83_074004.dat")
const δfuns = (δ00_0, δ11, δ20)
const dats  = (d00, d11, d20)
const ndata = d00.ndata + d11.ndata + d20.ndata
const nm = ["L$i" for i in 1:8]

function chi2_8(lec)
    s = 0.0
    for c in 1:3
        d = dats[c]; δ = δfuns[c]
        @inbounds for i in 1:d.ndata
            s += (sind(d.y[i] - δ(d.x[i], lec)) / (d.err[i]*π/180))^2
        end
    end
    return s
end

struct IAMPoint; chan::Int; x::Float64; y::Float64; err::Float64; end
const pts = IAMPoint[]
for c in 1:3, i in 1:dats[c].ndata
    push!(pts, IAMPoint(c, dats[c].x[i], dats[c].y[i], dats[c].err[i]))
end

# The `refit` contract: fit a resampled subset, warm-started from `start`, and
# return the fitted parameter VECTOR (NaNs for an invalid fit so the resampling
# dispatch drops it). Strategy-1 keeps discovery fast.
function iam_refit(subpts, start)
    function chi2r(lec)
        s = 0.0
        @inbounds for p in subpts
            s += (sind(p.y - δfuns[p.chan](p.x, lec)) / (p.err*π/180))^2
        end
        return s
    end
    fm = migrad(JuMinuit.CostFunction(chi2r, 1.0), start, fill(1e-6, 8);
                strategy = JuMinuit.Strategy(1))
    return JuMinuit.is_valid(fm) ? collect(fm.state.parameters.x) : fill(NaN, 8)
end

println("\n", "="^78)
println("find_deeper_minimum on the IAM ππ fit (data-resampling dispatch)")
println("="^78)
println("data points = $ndata,  free LECs = 8,  DOF = $(ndata - 8)")

# ── [1] cold single MIGRAD — lands in a shallow basin ────────────────────────
m = Minuit(chi2_8, collect(lecr0); names = nm, errors = fill(1e-6, 8), strategy = 1)
migrad!(m); hesse(m)
χ2_cold = m.fval
println("\n[1] cold Strategy-1 MIGRAD from the published LEC start:")
println("    χ² = ", round(χ2_cold; digits=3),
        "  (χ²/dof = ", round(χ2_cold/(ndata-8); digits=2), ")  valid = ", m.valid)

# ── [2] find_deeper_minimum — bootstrap-driven basin hopping ─────────────────
println("\n[2] find_deeper_minimum(m, iam_refit, pts; n_discovery=$NDISC, seed=$SEED):")
m_deep = find_deeper_minimum(m, iam_refit, pts;
                             n_discovery = NDISC, max_rounds = 6, seed = SEED, verbose = true)
χ2_deep = m_deep.fval

# ── [3] result ───────────────────────────────────────────────────────────────
println("\n[3] result:")
println("    χ² = ", round(χ2_deep; digits=3),
        "  (χ²/dof = ", round(χ2_deep/(ndata-8); digits=2), ")  valid = ", m_deep.valid)
println("    Δχ² over the cold fit = ", round(χ2_cold - χ2_deep; digits=3))
if χ2_deep < χ2_cold - 1
    println("    ✓ find_deeper_minimum escaped the shallow basin: χ² dropped ",
            round(χ2_cold; digits=1), " → ", round(χ2_deep; digits=1),
            " (Δ = ", round(χ2_cold - χ2_deep; digits=1), ").")
    println("    Now run error analysis (HESSE / MINOS / MC-Δχ²) at THIS minimum, not the cold one.")
else
    println("    (no deeper basin found from this start — raise n_discovery or try another seed.)")
end
println("="^78)
