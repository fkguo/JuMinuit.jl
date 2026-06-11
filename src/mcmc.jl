# SPDX-License-Identifier: LGPL-2.1-or-later

# ─────────────────────────────────────────────────────────────────────────────
# mcmc.jl — likelihood-ensemble MCMC (random-walk Metropolis) + quantile bands.
#
# The second leg of the error-analysis triangulation
#
#     profile extremization  ↔  likelihood-ensemble quantiles  ↔  MINOS
#
# (see docs/src/error_analysis.md): sample parameter sets from the
# likelihood L ∝ exp(−χ²/2) — the EXACT FCN, re-evaluated at every step,
# never a Gaussian surrogate — then turn any derived quantity f(θ) into
# marginal quantiles ([`quantiles`](@ref)) or pointwise quantile bands
# ([`quantile_band`](@ref)) over the ensemble.
#
# This is deliberately a *plain* single-chain random-walk Metropolis:
#
#   • The proposal (HESSE covariance / per-coordinate errors / explicit
#     user steps) affects MIXING EFFICIENCY only, never the stationary
#     distribution — any symmetric proposal targets exp(−Δf/(2·up)).
#     Contrast `get_contours_samples`, where an under-covering proposal
#     biases the reported region extents.
#   • Parameter limits are enforced by REJECTION, so the chain samples
#     the likelihood truncated to the allowed box. Posterior mass piling
#     up one-sidedly at an active boundary (a best fit sitting at g ≥ 0,
#     say) is the correct truncated marginal there — NOT a sampler bug,
#     and the reason a quantile band may legitimately exclude the best
#     fit (mode ≠ median).
#   • Step-size adaptation (`target_accept`) runs ONLY during burn-in;
#     the kept chain has a fixed kernel (exact detailed balance).
#
# iminuit has no native analogue (Python users bolt on emcee). The
# defaults encode validated field practice from real coupled-channel
# analyses: proposal ≈ 0.25–0.35 × HESSE σ → acceptance 0.2–0.3;
# 52 k steps, burn 2 k, thin 25 → 2000 kept sets, saved to disk
# (`save_ensemble`) and reused for every later derived quantity.
# ─────────────────────────────────────────────────────────────────────────────

"""
    LikelihoodEnsemble

Likelihood-weighted parameter ensemble produced by [`mcmc_sample`](@ref)
(or reloaded by [`load_ensemble`](@ref)). One row per kept Metropolis
step, in **full external** coordinates (fixed parameters appear as
constant columns), so any user model `f(θ_full)` evaluates directly on a
row.

# Fields

- `samples::Matrix{Float64}` — `nkept × npar` kept parameter sets (row =
  one set, columns ordered as `m.parameters`).
- `fvals::Vector{Float64}` — the FCN value (χ² or negative log-likelihood,
  whatever the fit minimized) at each kept set. Per-sample
  χ²-equivalent displacements are `(fvals .- fbest) ./ up`.
- `names::Vector{String}` — all parameter names (length `npar`).
- `free::Vector{Bool}` — which columns were varied by the chain.
- `best::Vector{Float64}` — the chain's start point (the fit's best
  values, full external vector).
- `fbest::Float64` — the FCN evaluated at `best` when the chain started.
- `up::Float64` — the fit's `errordef` (1 for χ², 0.5 for `−log L`).
- `acceptance::Float64` — accepted fraction of the post-burn-in proposals
  (healthy random-walk Metropolis: ≈ 0.2–0.4).
- `nsteps::Int`, `burn::Int`, `thin::Int` — the chain settings;
  `nkept = (nsteps − burn) ÷ thin`.
- `scale::Float64` — the final proposal scale (after any burn-in
  adaptation; equals the input `scale` when `target_accept` was not used).
- `proposal::Symbol` — the proposal actually used: `:hesse`, `:errors`,
  `:steps` (explicit per-coordinate σ) or `:matrix` (explicit covariance);
  `:unknown` for ensembles loaded from a foreign file.
- `seed::Union{Nothing,UInt64}` — the RNG seed, when one was given.

# Behaves like a collection of parameter vectors

`length(ens)` is the number of kept sets; `ens[i]` returns the `i`-th
parameter vector (a copy); iteration yields the rows, so
`[f(θ) for θ in ens]` evaluates a derived quantity over the whole
ensemble (which is exactly what [`quantiles`](@ref) and
[`quantile_band`](@ref) automate).

A chain can wander below the fit minimum (multi-modal landscape, or an
unconverged fit): `minimum(ens.fvals) < ens.fbest` is then a *finding*,
not an error — re-minimize (see `find_deeper_minimum`) before quoting
errors about the old minimum. The `show` method flags this.
"""
struct LikelihoodEnsemble
    samples::Matrix{Float64}
    fvals::Vector{Float64}
    names::Vector{String}
    free::Vector{Bool}
    best::Vector{Float64}
    fbest::Float64
    up::Float64
    acceptance::Float64
    nsteps::Int
    burn::Int
    thin::Int
    scale::Float64
    proposal::Symbol
    seed::Union{Nothing,UInt64}

    function LikelihoodEnsemble(samples, fvals, names, free, best, fbest, up,
                                acceptance, nsteps, burn, thin, scale, proposal,
                                seed)
        nk, np = size(samples)
        length(fvals) == nk ||
            throw(DimensionMismatch("fvals length $(length(fvals)) ≠ number of samples $nk"))
        length(names) == np ||
            throw(DimensionMismatch("names length $(length(names)) ≠ number of parameters $np"))
        length(free) == np ||
            throw(DimensionMismatch("free length $(length(free)) ≠ number of parameters $np"))
        length(best) == np ||
            throw(DimensionMismatch("best length $(length(best)) ≠ number of parameters $np"))
        return new(samples, fvals, names, free, best, fbest, up, acceptance,
                   nsteps, burn, thin, scale, proposal, seed)
    end
end

Base.length(e::LikelihoodEnsemble) = size(e.samples, 1)
Base.firstindex(e::LikelihoodEnsemble) = 1
Base.lastindex(e::LikelihoodEnsemble) = length(e)
Base.getindex(e::LikelihoodEnsemble, i::Integer) = e.samples[i, :]
Base.eltype(::Type{LikelihoodEnsemble}) = Vector{Float64}
function Base.iterate(e::LikelihoodEnsemble, i::Int = 1)
    i > length(e) && return nothing
    return e.samples[i, :], i + 1
end

function Base.show(io::IO, e::LikelihoodEnsemble)
    print(io, "LikelihoodEnsemble(", length(e), " samples, ",
          count(e.free), " free parameters)")
end

function Base.show(io::IO, ::MIME"text/plain", e::LikelihoodEnsemble)
    nk = length(e)
    nfree = count(e.free)
    ntot = length(e.names)
    print(io, "LikelihoodEnsemble: ", nk, " samples × ", nfree, " free parameter",
          nfree == 1 ? "" : "s")
    nfree < ntot && print(io, " (of ", ntot, ")")
    println(io)
    if isfinite(e.acceptance)
        @printf(io, "  acceptance %.3f   proposal :%s   scale %.4g   (nsteps %d, burn %d, thin %d)\n",
                e.acceptance, e.proposal, e.scale, e.nsteps, e.burn, e.thin)
    end
    if nk > 0 && isfinite(e.fbest) && isfinite(e.up) && e.up > 0
        Δ = (e.fvals .- e.fbest) ./ e.up
        q16, q50, q84 = (quantile(Δ, p) for p in (0.16, 0.5, 0.84))
        @printf(io, "  fval at start %.6g;  ensemble Δχ²-equivalent 16/50/84%%: %.3g / %.3g / %.3g\n",
                e.fbest, q16, q50, q84)
        lowest = minimum(e.fvals)
        if lowest < e.fbest - 1e-9 * max(1.0, abs(e.fbest))
            @printf(io, "  ⚠ chain found fval %.6g BELOW the start point %.6g — the fit may not be \
the global minimum (see `find_deeper_minimum`)\n", lowest, e.fbest)
        end
    end
end

# Same unreliability test as the `get_contours_samples` warning path: an
# invalid minimum, a forced-positive-definite covariance, or a failed /
# non-pos-def HESSE all mean Σ (and its diagonal) cannot be trusted as a
# *statement of errors* — though even then it usually remains a usable
# proposal SCALE, which is all the Metropolis chain needs.
function _covariance_unreliable(m::Minuit)
    fm = m.fmin.internal
    cov_status = fm.state.error.status
    return !is_valid(m.fmin) || fm.made_pos_def ||
           cov_status == MnHesseFailed || cov_status == MnMadePosDef ||
           cov_status == MnNotPosDef
end

"""
    mcmc_sample(m::Minuit; nsteps=52_000, burn=2_000, thin=25,
                proposal=:hesse, scale=0.3, target_accept=nothing,
                adapt_every=100, seed=nothing, rng=nothing, warn=true)
        -> LikelihoodEnsemble

Sample parameter sets from the fit's **likelihood** `L(θ) ∝
exp(−fcn(θ) / (2·up))` with a random-walk Metropolis chain on the **exact
FCN** — χ²: `exp(−Δχ²/2)`; negative log-likelihood (`up = 0.5`):
`exp(−Δ(−log L))` — started at the current best fit. The returned
[`LikelihoodEnsemble`](@ref) feeds [`quantiles`](@ref) /
[`quantile_band`](@ref) (marginal quantile intervals and pointwise bands
for any derived quantity) and can be stored with [`save_ensemble`](@ref)
as a reusable error set.

This is the **likelihood-ensemble leg** of the error-triangulation
(profile/MINOS extremization ↔ ensemble quantiles); iminuit has no native
analogue (Python users attach `emcee`). Call `migrad!(m)` (and ideally
`hesse!(m)`, for a well-shaped proposal) first.

# NOT the same as `get_contours_samples`

[`get_contours_samples`](@ref) samples the **confidence region** `Δχ² ≤
delta_chisq(cl, ndof)` — a hard cut, giving region *extents*. This
function runs a **posterior chain**: no Δχ² cut at all; samples
concentrate where the likelihood mass is, i.e. at Δχ² ≈ `n_free` (the
high-dimensional *volume effect* — for 9 free parameters
`P(Δχ² ≤ 1) ≈ 5.6e-4`, so a likelihood chain essentially never visits the
Δχ² ≤ 1 region, and does not need to). Use the region sampler for joint
confidence regions; use this chain for **likelihood-weighted quantiles of
derived quantities**.

# Marginal quantiles vs profile envelopes — read before quoting bands

A `(16%, 84%)` quantile band built from this ensemble is a **marginal
(posterior-mass) construction**. The profile construction (`mnprofile`,
MINOS, or constrained extremization of a derived quantity over
`Δχ² ≤ 1`) **contains the best fit by construction**; a marginal quantile
band need not. The two agree in the near-Gaussian interior and separate
legitimately at parameter limits: the truncation piles the posterior mass
one-sidedly, so the band can shift away from a best fit sitting on the
boundary (mode ≠ median). That separation is a *property of the
constructions*, not a sampler failure, and more samples will not make it
go away. See `docs/src/error_analysis.md` for the full comparison table.

# Algorithm

Plain Metropolis: propose `θ' = θ + scale·δ` with symmetric `δ`, accept
with probability `min(1, exp(−(fcn(θ′) − fcn(θ)) / (2·up)))`, record every
`thin`-th state after the first `burn` steps (`nkept = (nsteps − burn) ÷
thin`). Fixed parameters never move. Proposals **outside parameter
limits are rejected** before the FCN is called — the chain samples the
likelihood truncated to the allowed box (see above). Non-finite FCN
values are likewise rejected, never accepted. The chain bypasses the
fit's call counter (`m.nfcn` is untouched) and never mutates `m`.

The proposal shape only affects **mixing efficiency**, never what the
chain converges to (any symmetric proposal has the same stationary
distribution) — so an imperfect Σ is far less dangerous here than in
`get_contours_samples`, where proposal under-coverage biases the region.

# Keyword arguments

- `nsteps=52_000`, `burn=2_000`, `thin=25` — chain length, discarded
  burn-in, and keep-every-`thin` stride; the defaults yield 2000 kept
  sets (a field-validated recipe). Thinning mainly buys decorrelation;
  for quantile work ~1–2 k kept sets are usually plenty.
- `proposal=:hesse` — the proposal step shape, one of
  - `:hesse` — multivariate Gaussian with the fit covariance (correlations
    included). Falls back to `:errors` with a warning when no covariance
    is available or it looks unreliable (invalid fit / forced pos-def /
    failed HESSE);
  - `:errors` — independent per-coordinate Gaussians with σ = the fit's
    parabolic errors as frozen in the fit result (`m.fmin.ext_errors` —
    numerically `m.errors` right after the fit) — the classic hand-rolled
    choice;
  - an `AbstractVector` of length `n_free` — explicit per-free-parameter
    σ (external units). The escape hatch when both the covariance *and*
    the parabolic errors are meaningless, e.g. a parameter pinned at a
    limit;
  - an `AbstractMatrix` (`n_free × n_free`) — explicit proposal
    covariance.
- `scale=0.3` — overall step multiplier: step σ = `scale ×` (proposal σ).
  Field experience: `0.25–0.35 ×` HESSE σ gives acceptance ≈ 0.2–0.3 for
  ~10 free parameters. (Random-walk theory suggests `≈ 2.38/√n_free` for
  a perfectly matched `:hesse` proposal; in practice start at 0.3 and let
  `target_accept` tune. In LOW dimension the same `scale` accepts much
  more — ≈ 0.8 for a 3-parameter Gaussian fit — which is conservative,
  not a misconfiguration; thinning absorbs the extra correlation.)
- `target_accept=nothing` — when set (e.g. `0.25`), the scale is adapted
  toward that acceptance **during burn-in only** (every `adapt_every`
  steps, multiplicative update clamped to ×[0.5, 2] per round), then
  frozen, so the kept chain has a fixed kernel. Requires `burn ≥
  adapt_every`.
- `seed` — make the chain reproducible (`seed=11` etc.); `rng` — supply
  your own `AbstractRNG` instead (mutually exclusive with `seed`).
- `warn=true` — emit tuning warnings (unreliable-covariance fallback,
  post-burn acceptance < 0.05 or > 0.9).

# Returns

A [`LikelihoodEnsemble`](@ref): kept sets in `.samples` (full external
vectors, fixed columns constant), FCN values in `.fvals`, post-burn
`.acceptance`, the final `.scale`, and the chain metadata.

# Example

```julia
m = Minuit(chi2, x0; names = names, limits = limits)
migrad!(m); hesse!(m)

ens = mcmc_sample(m; seed = 11)                    # 52 k steps → 2000 sets
ens.acceptance                                     # aim for ≈ 0.2–0.3

# scalar derived quantity: 16/50/84% quantiles
q16, q50, q84 = quantiles(ens, θ -> θ[1] - θ[2])

# pointwise 16–84% band of a curve
band = quantile_band(ens, (x, θ) -> model(x, θ), xgrid)

save_ensemble("ensemble.dat", ens; comment = "error set B")  # reusable
```

See also [`quantiles`](@ref), [`quantile_band`](@ref),
[`save_ensemble`](@ref), [`load_ensemble`](@ref),
[`get_contours_samples`](@ref), and `docs/src/error_analysis.md`.
"""
function mcmc_sample(m::Minuit;
                     nsteps::Integer = 52_000,
                     burn::Integer = 2_000,
                     thin::Integer = 25,
                     proposal::Union{Symbol,AbstractVector{<:Real},AbstractMatrix{<:Real}} = :hesse,
                     scale::Real = 0.3,
                     target_accept::Union{Nothing,Real} = nothing,
                     adapt_every::Integer = 100,
                     seed::Union{Nothing,Integer} = nothing,
                     rng::Union{Nothing,Random.AbstractRNG} = nothing,
                     warn::Bool = true)
    m.fmin === nothing &&
        throw(ArgumentError("Call `migrad!(m)` before `mcmc_sample(m)` — the chain starts at the best fit"))
    nsteps >= 1 || throw(ArgumentError("nsteps must be ≥ 1"))
    0 <= burn < nsteps || throw(ArgumentError("need 0 ≤ burn < nsteps, got burn=$burn, nsteps=$nsteps"))
    thin >= 1 || throw(ArgumentError("thin must be ≥ 1"))
    nsteps - burn >= thin ||
        throw(ArgumentError("no samples would be kept: need nsteps − burn ≥ thin " *
                            "(got nsteps=$nsteps, burn=$burn, thin=$thin)"))
    scale > 0 || throw(ArgumentError("scale must be > 0"))
    adapt_every >= 1 || throw(ArgumentError("adapt_every must be ≥ 1"))
    if target_accept !== nothing
        0 < target_accept < 1 ||
            throw(ArgumentError("target_accept must be in (0, 1), got $target_accept"))
        burn >= adapt_every ||
            throw(ArgumentError("target_accept adaptation runs during burn-in: need burn ≥ adapt_every " *
                                "(got burn=$burn, adapt_every=$adapt_every)"))
    end
    seed !== nothing && rng !== nothing &&
        throw(ArgumentError("pass either `seed` or `rng`, not both"))
    # The seed is round-tripped through the ensemble metadata as a UInt64,
    # so reject an out-of-range seed up front rather than running the whole
    # chain and only failing at the `UInt64(seed)` return (negative OR larger
    # than typemax(UInt64) — MersenneTwister accepts both, the metadata does not).
    seed !== nothing && !(0 <= seed <= typemax(UInt64)) &&
        throw(ArgumentError("seed must be in [0, typemax(UInt64)] (stored as a UInt64 in the " *
                            "ensemble metadata), got $seed"))
    rng_use = seed !== nothing ? Random.MersenneTwister(seed) :
              rng !== nothing ? rng : Random.default_rng()

    ntot = n_pars(m.params)
    free_idx = [i for i in 1:ntot if !is_fixed(m.params.pars[i])]
    nfree = length(free_idx)
    nfree >= 1 || throw(ArgumentError("no free parameters to sample"))

    best_full = collect(Float64, m.values)
    up = Float64(m.fcn.up)
    up > 0 || throw(ArgumentError("errordef (up) must be > 0, got $up"))

    # Per-free-parameter limits, NaN ⇒ unbounded on that side (same
    # convention as get_contours_samples). Enforced by rejection BEFORE
    # the FCN is called.
    lo_free = [m.params.pars[i].lower for i in free_idx]
    hi_free = [m.params.pars[i].upper for i in free_idx]

    # ── Resolve the proposal into either a per-coordinate σ vector or a
    #    covariance square-root factor S (step = scale·S·z). ──────────────
    local steps_vec::Union{Nothing,Vector{Float64}} = nothing
    local Sfac::Union{Nothing,Matrix{Float64}} = nothing
    local prop_used::Symbol
    fallback_to_errors = false
    if proposal isa AbstractVector
        length(proposal) == nfree ||
            throw(DimensionMismatch("proposal σ vector has length $(length(proposal)), " *
                                    "expected n_free = $nfree"))
        steps_vec = Float64.(collect(proposal))
        all(s -> isfinite(s) && s > 0, steps_vec) ||
            throw(ArgumentError("explicit proposal σ must all be finite and > 0"))
        prop_used = :steps
    elseif proposal isa AbstractMatrix
        size(proposal) == (nfree, nfree) ||
            throw(DimensionMismatch("proposal covariance is $(size(proposal)), " *
                                    "expected ($nfree, $nfree)"))
        # An explicit covariance is a deliberate user act: a non-PD matrix
        # would silently freeze the clamped directions (zero proposal spread)
        # — fail loudly instead.
        isposdef(Symmetric(Matrix{Float64}(proposal))) ||
            throw(ArgumentError("explicit proposal covariance must be positive-definite " *
                                "(for a diagonal proposal pass a σ vector instead)"))
        Sfac = _mvnormal_factor(proposal)
        prop_used = :matrix
    elseif proposal === :hesse
        cov = free_covariance(m.fmin)
        if cov === nothing
            fallback_to_errors = true
            warn && @warn "mcmc_sample: no covariance available for proposal=:hesse; " *
                          "falling back to the per-coordinate :errors proposal. " *
                          "Run `hesse!(m)` for a correlation-aware proposal."
        elseif _covariance_unreliable(m)
            fallback_to_errors = true
            warn && @warn """mcmc_sample: the fit covariance looks unreliable \
(is_valid=$(is_valid(m.fmin)), made_pos_def=$(m.fmin.internal.made_pos_def), \
cov_status=$(m.fmin.internal.state.error.status)); falling back to the per-coordinate \
:errors proposal. This only affects mixing efficiency, not what the chain converges to — \
but consider an explicit per-parameter σ vector (`proposal = [σ₁, σ₂, …]`) if the \
parabolic errors are also meaningless (e.g. a parameter at a limit)."""
        else
            Sfac = _mvnormal_factor(cov)
            prop_used = :hesse
        end
    elseif proposal === :errors
        fallback_to_errors = true   # not a fallback, the explicit request — same setup path
    else
        throw(ArgumentError("proposal must be :hesse, :errors, a σ vector, or a covariance matrix; got :$proposal"))
    end
    if fallback_to_errors
        errs = collect(Float64, m.fmin.ext_errors)[free_idx]
        all(e -> isfinite(e) && e > 0, errs) ||
            throw(ArgumentError("the fit's parabolic errors are not all finite and positive — " *
                                "pass an explicit σ vector: proposal = [σ₁, σ₂, …]"))
        steps_vec = errs
        prop_used = :errors
    end

    # χ²/FCN over the free coordinates: splice into a fresh full external
    # vector (the user FCN never sees our working buffers) and call the
    # raw user function — NOT the counting wrapper, so m.nfcn is untouched.
    userf = m.fcn.f
    eval_free = let base = best_full, fi = free_idx, f = userf
        q -> begin
            full = copy(base)
            @inbounds for (j, i) in enumerate(fi)
                full[i] = q[j]
            end
            return Float64(f(full))
        end
    end

    p = best_full[free_idx]                 # current chain position (free coords)
    c0 = eval_free(p)
    isfinite(c0) ||
        throw(ArgumentError("the FCN is not finite at the best-fit start point (fcn = $c0)"))
    fbest = c0

    nkept = (nsteps - burn) ÷ thin
    kept = Matrix{Float64}(undef, nkept, ntot)
    fvals = Vector{Float64}(undef, nkept)

    q = Vector{Float64}(undef, nfree)       # proposal buffer
    z = Vector{Float64}(undef, nfree)       # N(0,1) draws for :hesse/:matrix
    scale_cur = Float64(scale)
    ikeep = 0
    nacc_post = 0                            # accepted after burn-in
    npost = nsteps - burn
    block_acc = 0                            # burn-in adaptation bookkeeping
    block_n = 0
    adapting = target_accept !== nothing

    @inbounds for it in 1:nsteps
        # Propose q = p + scale·δ with symmetric δ.
        if steps_vec !== nothing
            for k in 1:nfree
                q[k] = p[k] + scale_cur * steps_vec[k] * randn(rng_use)
            end
        else
            for k in 1:nfree
                z[k] = randn(rng_use)
            end
            for k in 1:nfree
                acc = 0.0
                for l in 1:nfree
                    acc += Sfac[k, l] * z[l]
                end
                q[k] = p[k] + scale_cur * acc
            end
        end

        # Reject outside the limits BEFORE calling the FCN (truncated
        # target); otherwise standard Metropolis on exp(−f/(2up)).
        accepted = false
        inbox = true
        for k in 1:nfree
            (!isnan(lo_free[k]) && q[k] < lo_free[k]) && (inbox = false; break)
            (!isnan(hi_free[k]) && q[k] > hi_free[k]) && (inbox = false; break)
        end
        if inbox
            c1 = eval_free(q)
            if isfinite(c1) && (c1 < c0 || rand(rng_use) < exp(-(c1 - c0) / (2.0 * up)))
                copyto!(p, q)
                c0 = c1
                accepted = true
            end
        end

        if it > burn
            accepted && (nacc_post += 1)
            r = it - burn
            if r % thin == 0
                ikeep += 1
                for i in 1:ntot
                    kept[ikeep, i] = best_full[i]
                end
                for (j, i) in enumerate(free_idx)
                    kept[ikeep, i] = p[j]
                end
                fvals[ikeep] = c0
            end
        elseif adapting
            block_n += 1
            accepted && (block_acc += 1)
            if block_n == adapt_every
                rate = block_acc / adapt_every
                scale_cur *= clamp(rate / Float64(target_accept), 0.5, 2.0)
                scale_cur = clamp(scale_cur, 1e-10, 1e10)
                block_n = 0
                block_acc = 0
            end
        end
    end
    @assert ikeep == nkept

    acceptance = nacc_post / npost
    if warn
        acceptance < 0.05 &&
            @warn "mcmc_sample: post-burn acceptance $(round(acceptance; digits=3)) < 0.05 — " *
                  "the steps are too large (or the start point is poor); lower `scale` or " *
                  "set `target_accept = 0.25`."
        acceptance > 0.9 &&
            @warn "mcmc_sample: post-burn acceptance $(round(acceptance; digits=3)) > 0.9 — " *
                  "the steps are much too small, so consecutive samples are highly correlated; " *
                  "raise `scale` or set `target_accept = 0.25`."
    end

    names = [p_.name for p_ in m.params.pars]
    free = [!is_fixed(p_) for p_ in m.params.pars]
    return LikelihoodEnsemble(kept, fvals, names, free, best_full, fbest, up,
                              acceptance, Int(nsteps), Int(burn), Int(thin),
                              scale_cur, prop_used,
                              seed === nothing ? nothing : UInt64(seed))
end

# ─────────────────────────────────────────────────────────────────────────────
# Quantiles of derived quantities over the ensemble.
# ─────────────────────────────────────────────────────────────────────────────

# Quantiles of `vals` at probabilities `ps`, dropping non-finite entries
# (returns the number dropped). `vals` is sorted in place.
function _finite_quantiles!(vals::Vector{Float64}, ps)
    ndrop = count(!isfinite, vals)
    if ndrop > 0
        filter!(isfinite, vals)
        isempty(vals) &&
            throw(ArgumentError("the derived quantity is non-finite on every ensemble member"))
    end
    sort!(vals)
    return [quantile(vals, Float64(pp); sorted = true) for pp in ps], ndrop
end

_warn_dropped(fname::AbstractString, ndrop::Int, ntotal::Int) =
    ndrop > 0 && @warn "$fname: dropped $ndrop of $ntotal non-finite values of the derived quantity"

"""
    quantiles(ens::LikelihoodEnsemble, f; p=(0.16, 0.5, 0.84), warn=true)
        -> Vector{Float64}

Quantiles of a scalar derived quantity `f(θ)` over the likelihood
ensemble: evaluates `f` on every kept parameter set (full external
vector) and returns the quantiles at probabilities `p`, in order.

The default `p` gives the 16% / median / 84% triplet — a marginal
"1σ-equivalent" interval. Include `0` / `1` in `p` for the ensemble
minimum / maximum. Non-finite `f` values are dropped (with a warning).

These are **marginal (posterior-mass) quantiles** — see the
[`mcmc_sample`](@ref) docstring for how they relate to (and legitimately
differ from) profile/MINOS intervals, especially at parameter limits.

```julia
q16, q50, q84 = quantiles(ens, θ -> θ[2] - θ[1])
lo, hi        = quantiles(ens, θ -> θ[1]; p = (0.16, 0.84))
```

See also [`quantile_band`](@ref) for pointwise bands of a curve.
"""
function quantiles(ens::LikelihoodEnsemble, f; p = (0.16, 0.5, 0.84), warn::Bool = true)
    n = length(ens)
    n >= 1 || throw(ArgumentError("the ensemble is empty"))
    vals = Vector{Float64}(undef, n)
    # `f` gets a fresh COPY of each row (like `ens[i]`/iteration), so a
    # user model that mutates its argument cannot corrupt the ensemble.
    @inbounds for j in 1:n
        vals[j] = Float64(f(ens.samples[j, :]))
    end
    out, ndrop = _finite_quantiles!(vals, p)
    warn && _warn_dropped("quantiles", ndrop, n)
    return out
end

"""
    quantile_band(ens::LikelihoodEnsemble, f, xs; p=(0.16, 0.84),
                  curve=false, warn=true) -> Matrix{Float64}

Pointwise quantile band of a curve over the likelihood ensemble: at each
grid point `xs[i]`, the quantiles (at probabilities `p`) of
`f(xs[i], θ)` across all kept parameter sets θ. Returns a
`length(xs) × length(p)` matrix — with the default `p`, column 1 is the
16% (lower) and column 2 the 84% (upper) band edge.

- `f` — `f(x, θ)::Real` with `θ` the full external parameter vector.
  With `curve = true`, instead `f(θ)::AbstractVector` returning the whole
  curve over `xs` in one call (one `f` call per ensemble member — use
  this when evaluating the model pointwise repeats expensive shared
  work).
- `p` — quantile probabilities; e.g. `p = (0.025, 0.16, 0.5, 0.84, 0.975)`
  for median + 1σ + 2σ-equivalent bands.
- Non-finite values are dropped per grid point (with a warning).

This is a **marginal quantile band** (likelihood-mass construction), the
companion of — not a substitute for — the profile envelope band
(pointwise extremization of `f` over `Δχ² ≤ delta_chisq(cl, ndof)`): the
profile band contains the best-fit curve by construction, a quantile
band need not (boundary effects shift the mass one-sidedly; mode ≠
median). Quote which construction you used. See [`mcmc_sample`](@ref)
and `docs/src/error_analysis.md`.

```julia
band = quantile_band(ens, (x, θ) -> model(x, θ), xgrid)
plot(xgrid, mid; ribbon = (mid .- band[:, 1], band[:, 2] .- mid))

# expensive model: one call per ensemble member returns the whole curve
band = quantile_band(ens, θ -> model_curve(xgrid, θ), xgrid; curve = true)
```
"""
function quantile_band(ens::LikelihoodEnsemble, f, xs; p = (0.16, 0.84),
                       curve::Bool = false, warn::Bool = true)
    n = length(ens)
    n >= 1 || throw(ArgumentError("the ensemble is empty"))
    nx = length(xs)
    nx >= 1 || throw(ArgumentError("xs is empty"))
    np = length(p)
    Q = Matrix{Float64}(undef, nx, np)
    ndrop_tot = 0
    vals = Vector{Float64}(undef, n)
    if curve
        # One f call per ensemble member; f returns the curve over xs.
        # Row COPIES, not views — a mutating `f` must not corrupt the ensemble.
        Y = Matrix{Float64}(undef, n, nx)
        @inbounds for j in 1:n
            y = f(ens.samples[j, :])
            length(y) == nx ||
                throw(DimensionMismatch("curve=true: f(θ) returned length $(length(y)), expected length(xs) = $nx"))
            for i in 1:nx
                Y[j, i] = Float64(y[i])
            end
        end
        @inbounds for i in 1:nx
            for j in 1:n
                vals[j] = Y[j, i]
            end
            qs, ndrop = _finite_quantiles!(vals, p)
            ndrop_tot += ndrop
            for k in 1:np
                Q[i, k] = qs[k]
            end
            resize!(vals, n)
        end
    else
        # A FRESH row copy per call (like `quantiles`): a mutating `f` must
        # corrupt neither the stored ensemble NOR a later grid point's view
        # of the same member — caching one row per member and reusing it
        # across `xs` would leak in-place mutations from one grid point into
        # the next.
        @inbounds for (i, x) in enumerate(xs)
            for j in 1:n
                vals[j] = Float64(f(x, ens.samples[j, :]))
            end
            qs, ndrop = _finite_quantiles!(vals, p)
            ndrop_tot += ndrop
            for k in 1:np
                Q[i, k] = qs[k]
            end
            resize!(vals, n)
        end
    end
    warn && _warn_dropped("quantile_band", ndrop_tot, n * nx)
    return Q
end

# ─────────────────────────────────────────────────────────────────────────────
# Plain-text ensemble persistence (reusable error sets).
# ─────────────────────────────────────────────────────────────────────────────

"""
    save_ensemble(path_or_io, ens::LikelihoodEnsemble; comment="")

Write the ensemble as plain text: `#`-comment header (metadata +
optional user `comment`), then one line per kept set — the FCN value
followed by all parameter values, space-separated. Floats are written in
shortest round-trip form, so [`load_ensemble`](@ref) reproduces them
**exactly**.

The format is deliberately the time-honoured hand-rolled one (`# header`
+ `fval p₁ p₂ …` rows): existing files of that shape load fine, and the
saved file is directly consumable by gnuplot / numpy / a 5-line Julia
loop. Ensembles are worth saving: any *future* derived quantity gets its
band by evaluating over the stored sets — no re-sampling.

Because the columns are whitespace-separated, parameter **names must not
contain whitespace** (they are written to the `# names:` header); a name
with a space throws rather than corrupting the round-trip.

```julia
save_ensemble("ensemble_B.dat", ens; comment = "error set B, ρ = -0.295")
```
"""
function save_ensemble(io::IO, ens::LikelihoodEnsemble; comment::AbstractString = "")
    nk, np = size(ens.samples)
    # The format is whitespace-separated, so a parameter name containing
    # whitespace would not round-trip (load_ensemble would split it into
    # several names). Fail loudly rather than silently mangle the header.
    any(nm -> occursin(r"\s", nm), ens.names) &&
        throw(ArgumentError("save_ensemble: parameter names must not contain whitespace " *
                            "(the plain-text format is whitespace-separated); got $(ens.names)"))
    # Metadata first, user comment after: load_ensemble takes the FIRST
    # occurrence of each `# key: value`, so a comment line that happens to
    # look like one (e.g. "# up: 5 was too big") can never shadow the real
    # metadata on round-trip.
    println(io, "# JuMinuit LikelihoodEnsemble v1")
    println(io, "# names: ", join(ens.names, ' '))
    println(io, "# free: ", join(Int.(ens.free), ' '))
    println(io, "# best: ", join(ens.best, ' '))
    println(io, "# fbest: ", ens.fbest)
    println(io, "# up: ", ens.up)
    println(io, "# acceptance: ", ens.acceptance)
    println(io, "# nsteps: ", ens.nsteps)
    println(io, "# burn: ", ens.burn)
    println(io, "# thin: ", ens.thin)
    println(io, "# scale: ", ens.scale)
    println(io, "# proposal: ", ens.proposal)
    ens.seed === nothing || println(io, "# seed: ", ens.seed)
    for line in split(comment, '\n')
        isempty(line) || println(io, "# ", line)
    end
    println(io, "# cols: fval ", join(ens.names, ' '))
    for j in 1:nk
        print(io, ens.fvals[j])
        for i in 1:np
            print(io, ' ', ens.samples[j, i])
        end
        println(io)
    end
    return nothing
end

function save_ensemble(path::AbstractString, ens::LikelihoodEnsemble;
                       comment::AbstractString = "")
    open(io -> save_ensemble(io, ens; comment), path, "w")
    return nothing
end

"""
    load_ensemble(path_or_io; names=nothing, up=nothing) -> LikelihoodEnsemble

Read an ensemble written by [`save_ensemble`](@ref) — or any plain-text
file of the same shape: `#` comment lines, then one row per sample,
`fval p₁ p₂ …` space-separated (the classic hand-rolled ensemble
format).

Metadata found in `# key: value` header lines (as written by
`save_ensemble`) is restored; anything missing gets a placeholder
(`names = ["p1", …]`, `free` all `true`, `best`/`fbest`/`up`/… `NaN`,
`proposal = :unknown`). Override with the `names` / `up` keywords when
loading a foreign file. Quantile analysis ([`quantiles`](@ref) /
[`quantile_band`](@ref)) only needs the samples, so it works on foreign
files as-is.
"""
function load_ensemble(io::IO; names::Union{Nothing,AbstractVector} = nothing,
                       up::Union{Nothing,Real} = nothing)
    meta = Dict{String,String}()
    rows = Vector{Vector{Float64}}()
    ncol = -1
    for (lineno, raw) in enumerate(eachline(io))
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "#")
            mm = match(r"^#\s*([a-z_]+):\s*(.*)$", line)
            mm !== nothing && !haskey(meta, mm.captures[1]) &&
                (meta[mm.captures[1]] = String(strip(mm.captures[2])))
            continue
        end
        toks = split(line)
        vals = try
            [parse(Float64, t) for t in toks]
        catch
            throw(ArgumentError("load_ensemble: line $lineno is not a numeric row: $(repr(raw))"))
        end
        length(vals) >= 2 ||
            throw(ArgumentError("load_ensemble: line $lineno has $(length(vals)) column(s); " *
                                "need at least `fval p₁`"))
        ncol == -1 && (ncol = length(vals))
        length(vals) == ncol ||
            throw(ArgumentError("load_ensemble: line $lineno has $(length(vals)) columns, " *
                                "previous rows had $ncol"))
        push!(rows, vals)
    end
    isempty(rows) && throw(ArgumentError("load_ensemble: no data rows found"))
    nk = length(rows)
    np = ncol - 1
    samples = Matrix{Float64}(undef, nk, np)
    fvals = Vector{Float64}(undef, nk)
    @inbounds for j in 1:nk
        fvals[j] = rows[j][1]
        for i in 1:np
            samples[j, i] = rows[j][i + 1]
        end
    end

    _split_meta(key) = haskey(meta, key) ? split(meta[key]) : nothing
    _parsef(key, default) = haskey(meta, key) ? something(tryparse(Float64, meta[key]), default) : default
    _parsei(key, default) = haskey(meta, key) ? something(tryparse(Int, meta[key]), default) : default

    names_use = if names !== nothing
        length(names) == np ||
            throw(DimensionMismatch("names has length $(length(names)), file has $np parameter columns"))
        String.(collect(names))
    else
        toks = _split_meta("names")
        toks !== nothing && length(toks) == np ? String.(toks) : ["p$i" for i in 1:np]
    end
    free = let toks = _split_meta("free")
        if toks !== nothing && length(toks) == np
            [t != "0" for t in toks]
        else
            fill(true, np)
        end
    end
    best = let toks = _split_meta("best")
        if toks !== nothing && length(toks) == np
            [something(tryparse(Float64, t), NaN) for t in toks]
        else
            fill(NaN, np)
        end
    end
    up_use = up !== nothing ? Float64(up) : _parsef("up", NaN)
    proposal = haskey(meta, "proposal") ? Symbol(meta["proposal"]) : :unknown
    seed = haskey(meta, "seed") ? tryparse(UInt64, meta["seed"]) : nothing

    return LikelihoodEnsemble(samples, fvals, names_use, free, best,
                              _parsef("fbest", NaN), up_use,
                              _parsef("acceptance", NaN),
                              _parsei("nsteps", 0), _parsei("burn", 0),
                              _parsei("thin", 0), _parsef("scale", NaN),
                              proposal, seed)
end

function load_ensemble(path::AbstractString; kwargs...)
    return open(io -> load_ensemble(io; kwargs...), path, "r")
end
