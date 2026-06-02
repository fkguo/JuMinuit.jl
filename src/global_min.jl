# SPDX-License-Identifier: LGPL-2.1-or-later
#
# Basin-hopping global minimisation for multi-basin objectives.
#
# MIGRAD (like every local optimiser) converges to whatever basin its start
# point drains into. On an ill-conditioned, multi-basin surface — IAM ππ being
# the worked example in `BenchmarkExamples/IAM_2Pformfactor` — that basin is
# often NOT the global one, and any error analysis done there is meaningless.
# `find_global_minimum` automates the "restart, find a deeper basin, adopt it,
# repeat" loop a user would otherwise run by hand.
#
# This is the global-SEARCH counterpart to `find_solution_modes` (which CLUSTERS
# an already-sampled set into distinct solutions). Use this to find the true
# minimum first, then the usual error analysis (HESSE / MINOS / get_contours_
# samples) at the minimum it returns.

"""
    find_global_minimum(fcn, x0, errors; kwargs...) -> FunctionMinimum

Basin-hopping global search for a multi-basin objective. Starting from a MIGRAD
fit at `x0`, repeatedly draw `n_restarts` perturbed restarts around the current
best (each coordinate jittered by `perturb · scaleᵢ · randn`, with
`scaleᵢ = max(|xᵢ|, |errorᵢ|, abs_floor)`), MIGRAD each, and **adopt any deeper
valid minimum**. A round that finds no improvement means the search has
converged; otherwise it stops after `max_rounds`. Returns the deepest
[`FunctionMinimum`](@ref) found.

`fcn` may be a plain callable (wrapped in `CostFunction(fcn, up)`) or an
[`AbstractCostFunction`](@ref). `x0`/`errors` are the usual MIGRAD start point
and step sizes.

!!! warning "Unbounded only — and check validity"
    This routine fits through the **unbounded** MIGRAD path and **ignores
    parameter limits**: fold any bounds into `fcn` (a penalty, or a smooth
    reparameterisation) before calling. The returned `FunctionMinimum` can be
    invalid (e.g. if every restart failed) — always check `is_valid(result)`
    before using it.

# Keyword arguments

- `n_restarts::Integer = 24` — perturbed restarts per round (must be ≥ 1).
- `perturb::Real = 1.0` — exploration radius, as a multiple of each parameter's
  scale. Larger ⇒ jumps farther (more likely to escape a basin, at one full
  re-fit per restart). **The key knob to tune on a hard surface.**
- `abs_floor::Real = 0.0` — an absolute lower bound on each coordinate's jitter
  scale. Raise it if a parameter sits near 0 with a tiny step (otherwise
  `scaleᵢ → 0` and that coordinate is never explored).
- `max_rounds::Integer = 6` — stop after this many improvement rounds (≥ 1).
- `strategy = Strategy(1)` — MIGRAD strategy for every (re-)fit.
- `maxfcn::Union{Integer,Nothing} = nothing` — per-fit MIGRAD call budget
  (`nothing` ⇒ MIGRAD's default `200 + 100n + 5n²`). Raise it for expensive FCNs
  whose restarts need many calls to converge.
- `min_improvement::Real = 1e-3` — a restart must beat the current best χ² by
  more than this to be adopted (guards against same-basin numerical jitter).
  (Note: this is the *adoption margin*, NOT MIGRAD's EDM `tol`.)
- `up::Real = 1.0` — error definition, when `fcn` is a bare callable.
- `seed::Union{Integer,Nothing} = nothing` — RNG seed for reproducible restarts.
- `verbose::Bool = false` — log the best χ² per round.

!!! note "Not a global-optimum guarantee"
    Basin-hopping is a heuristic: it finds a *deeper* basin when its restarts
    land in one, but cannot prove the result is global. Raise `n_restarts` /
    `perturb` / `max_rounds` for a more thorough search, and cross-check by
    re-running from independent seeds.

# Example

```julia
fm = find_global_minimum(chi2, x0, errs; n_restarts = 40, perturb = 1.5, seed = 1)
is_valid(fm) || error("global search failed")
m = Minuit(chi2, values(fm); names = pnames)   # error analysis at the minimum
migrad!(m); hesse(m)
```

See also [`find_solution_modes`](@ref) (cluster sampled solutions into modes).
"""
function find_global_minimum(cf::AbstractCostFunction, x0::AbstractVector, errors::AbstractVector;
        n_restarts::Integer = 24, perturb::Real = 1.0, abs_floor::Real = 0.0,
        max_rounds::Integer = 6, strategy = Strategy(1),
        maxfcn::Union{Integer,Nothing} = nothing, min_improvement::Real = 1e-3,
        seed::Union{Integer,Nothing} = nothing, verbose::Bool = false)
    n_restarts >= 1 || throw(ArgumentError("find_global_minimum: n_restarts must be ≥ 1"))
    max_rounds >= 1 || throw(ArgumentError("find_global_minimum: max_rounds must be ≥ 1"))
    perturb > 0 || throw(ArgumentError("find_global_minimum: perturb must be > 0"))
    min_improvement >= 0 || throw(ArgumentError("find_global_minimum: min_improvement must be ≥ 0"))
    abs_floor >= 0 || throw(ArgumentError("find_global_minimum: abs_floor must be ≥ 0"))
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    errs = collect(Float64, errors)

    # Each fit gets a fresh call budget: `migrad` compares the cost function's
    # CUMULATIVE `nfcn` to `maxfcn`, so without resetting, later restarts on the
    # shared `cf` would start already over the limit and bail immediately.
    _fit(x) = (reset_ncalls!(cf); migrad(cf, x, errs; strategy = strategy, maxfcn = maxfcn))

    best = _fit(collect(Float64, x0))
    for round in 1:max_rounds
        bx = collect(Float64, best.state.parameters.x)
        scale = [max(abs(bx[i]), abs(errs[i]), abs_floor, eps()) for i in eachindex(bx)]
        improved = false
        for _ in 1:n_restarts
            x = bx .+ perturb .* scale .* randn(rng, length(bx))
            # A wild jitter can push a constrained FCN into a throwing region
            # (log of a negative, a singular matrix, …); skip that restart
            # rather than aborting the whole search.
            fm = try
                _fit(x)
            catch err
                err isa Union{DomainError,BoundsError,SingularException,ArgumentError,DivideError} || rethrow()
                continue
            end
            # NB: restarts in a round all jitter around this round's `bx` (the
            # current best is re-centred only at the next round boundary) — the
            # standard basin-hopping choice.
            if is_valid(fm) && isfinite(fval(fm)) && fval(fm) < fval(best) - min_improvement
                best = fm
                improved = true
            end
        end
        verbose && @info "find_global_minimum" round χ² = fval(best) improved
        improved || break
    end
    return best
end

find_global_minimum(f, x0::AbstractVector, errors::AbstractVector; up::Real = 1.0, kwargs...) =
    find_global_minimum(CostFunction(f, up), x0, errors; kwargs...)
