# SPDX-License-Identifier: LGPL-2.1-or-later
using JuMinuit
using Test
using Logging

@testset "find_deeper_minimum (basin-hopping local-escape search)" begin
    # Double well in x[1] (minima ≈ ±1), tilted by +0.4·x[1] so the x[1]≈−1 well
    # is DEEPER (f ≈ −0.41) than the x[1]≈+1 well (f ≈ +0.39); x[2] is a simple
    # quadratic. A plain MIGRAD started on the +1 side stays in the shallow basin;
    # find_deeper_minimum must escape to the deep one. (A wide search —
    # n_restarts=80, perturb=2 — makes the escape robust, not seed-lucky.)
    f(x) = (x[1]^2 - 1)^2 + 0.4 * x[1] + x[2]^2
    shallow = 0.4    # ≈ f at the +1 well

    fm = find_deeper_minimum(f, [1.0, 0.5], [0.3, 0.3];
                             n_restarts = 80, perturb = 2.0, seed = 1)
    @test JuMinuit.is_valid(fm)
    @test fm.state.parameters.x[1] < 0          # escaped to the deeper (−1) well
    @test JuMinuit.fval(fm) < shallow - 0.4     # and it is genuinely deeper
    @test values(fm) ≈ collect(fm.state.parameters.x)   # public accessor agrees

    # reproducible (same seed ⇒ same result)
    fm2 = find_deeper_minimum(f, [1.0, 0.5], [0.3, 0.3];
                              n_restarts = 80, perturb = 2.0, seed = 1)
    @test fm2.state.parameters.x ≈ fm.state.parameters.x

    # the AbstractCostFunction form agrees with the bare-callable form
    fmc = find_deeper_minimum(CostFunction(f, 1.0), [1.0, 0.5], [0.3, 0.3];
                              n_restarts = 80, perturb = 2.0, seed = 1)
    @test fmc.state.parameters.x ≈ fm.state.parameters.x

    # A throwing FCN must NOT abort the search: log(x[1]) throws for x[1] ≤ 0,
    # and the wide jitter will probe there — those restarts are skipped.
    g(x) = (log(x[1]))^2 + (x[1] - 2.0)^2 + x[2]^2
    fmg = find_deeper_minimum(g, [2.0, 0.0], [0.5, 0.5];
                              n_restarts = 40, perturb = 1.0, seed = 3)
    @test JuMinuit.is_valid(fmg)
    @test fmg.state.parameters.x[1] > 0         # converged in the valid domain

    # argument validation
    @test_throws ArgumentError find_deeper_minimum(f, [0.0], [0.1]; n_restarts = 0)
    @test_throws ArgumentError find_deeper_minimum(f, [0.0], [0.1]; perturb = 0.0)
    @test_throws ArgumentError find_deeper_minimum(f, [0.0], [0.1]; max_rounds = 0)
    @test_throws ArgumentError find_deeper_minimum(f, [0.0], [0.1]; min_improvement = -1.0)

    # the deprecated v0.3.1 name still forwards to find_deeper_minimum
    # (depwarn visibility depends on the --depwarn flag, so just check forwarding)
    fmd = find_global_minimum(f, [1.0, 0.5], [0.3, 0.3]; n_restarts = 80, perturb = 2.0, seed = 1)
    @test fmd.state.parameters.x ≈ fm.state.parameters.x
end

# ─────────────────────────────────────────────────────────────────────────────
# New dispatches added in v0.3.2
# ─────────────────────────────────────────────────────────────────────────────

@testset "find_deeper_minimum — Minuit convenience dispatch (perturbation)" begin
    # Passing a converged Minuit should give the same result as (cf, x0, errors).
    f(x) = (x[1]^2 - 1)^2 + 0.4 * x[1] + x[2]^2
    cf   = CostFunction(f, 1.0)
    x0   = [1.0, 0.5]; errs = [0.3, 0.3]

    m_ref = Minuit(f, x0; errors = errs, strategy = 1)
    migrad!(m_ref)

    # Both paths: same seed → same trajectory → same minimum.
    fm_via_m  = find_deeper_minimum(m_ref;    n_restarts = 80, perturb = 2.0, seed = 7)
    fm_direct = find_deeper_minimum(cf, x0, errs; n_restarts = 80, perturb = 2.0, seed = 7)

    @test fm_via_m  isa JuMinuit.FunctionMinimum
    @test fm_direct isa JuMinuit.FunctionMinimum
    @test JuMinuit.is_valid(fm_via_m)
    # Same basin and very close fval (tiny numerical drift from one extra initial migrad call).
    @test JuMinuit.fval(fm_via_m) ≈ JuMinuit.fval(fm_direct) atol = 0.01
    @test values(fm_via_m)[1] * values(fm_direct)[1] > 0  # same sign ⇒ same basin
end

@testset "find_deeper_minimum — resampling dispatches" begin
    # ── Single-basin fixture ──────────────────────────────────────────────────
    # chi2(p) = Σ_i (y_i - p[1])^2 / 0.01, data all at y=2.
    # Every bootstrap resample converges to p≈2 → suitability check fires.
    pts_sb = fill(2.0, 40)

    function chi2_sb(p, d = pts_sb)
        s = 0.0
        for y in d; s += (y - p[1])^2 / 0.01; end
        return s
    end
    refit_sb = (subdata, start) -> begin
        fm = migrad(CostFunction(p -> chi2_sb(p, subdata), 1.0),
                    start, [0.1]; strategy = Strategy(1))
        JuMinuit.is_valid(fm) ? collect(Float64, values(fm)) : fill(NaN, length(start))
    end

    m_sb = Minuit(chi2_sb, [1.9]; errors = [0.1])
    migrad!(m_sb); hesse(m_sb)

    # ── pre-fitted dispatch: suitability warning + return m unchanged ─────────
    @testset "suitability check — single-basin warns and returns m" begin
        m_out = @test_logs((:warn, r"No deeper basin|no deeper basin"),
                           min_level = Logging.Warn,
                           find_deeper_minimum(m_sb, refit_sb, pts_sb;
                                               n_discovery = 10, seed = 1))
        @test m_out isa Minuit
        # Returned the same object (no adoption happened).
        @test m_out === m_sb
    end

    # ── Argument validation ───────────────────────────────────────────────────
    @testset "argument validation — resampling dispatch" begin
        @test_throws ArgumentError find_deeper_minimum(m_sb, refit_sb, pts_sb;
                                                       n_discovery = 1)
        @test_throws ArgumentError find_deeper_minimum(m_sb, refit_sb, pts_sb;
                                                       max_rounds = 0)
        @test_throws ArgumentError find_deeper_minimum(m_sb, refit_sb, pts_sb;
                                                       min_improvement = -0.1)
        @test_throws ArgumentError find_deeper_minimum(m_sb, refit_sb, [2.0];  # length < 2
                                                       n_discovery = 2)
    end

    # ── fresh-start delegates to pre-fitted: same result ─────────────────────
    @testset "fresh-start delegates to pre-fitted (single-basin, warns once)" begin
        cf_sb  = CostFunction(chi2_sb, 1.0)
        # Both should fire suitability warning and return a Minuit near chi2≈0.
        m_pre   = @test_logs((:warn, r"No deeper basin|no deeper basin"),
                             min_level = Logging.Warn,
                             find_deeper_minimum(m_sb, refit_sb, pts_sb;
                                                 n_discovery = 10, seed = 2))
        m_fresh = @test_logs((:warn, r"No deeper basin|no deeper basin"),
                             min_level = Logging.Warn,
                             find_deeper_minimum(cf_sb, [1.9], [0.1], refit_sb, pts_sb;
                                                 n_discovery = 10, seed = 2))
        @test m_pre   isa Minuit
        @test m_fresh isa Minuit
        @test m_pre.fval  ≈ m_fresh.fval  atol = 1e-4
    end

    # ── plain-callable wrapper reaches the same path ──────────────────────────
    @testset "plain-callable wrapper for resampling" begin
        m_plain = @test_logs((:warn, r"No deeper basin|no deeper basin"),
                             min_level = Logging.Warn,
                             find_deeper_minimum(chi2_sb, [1.9], [0.1], refit_sb, pts_sb;
                                                 n_discovery = 10, seed = 3))
        @test m_plain isa Minuit
    end

    # ── disambiguator: (Minuit, AbstractVector, AbstractVector) → ArgumentError ─
    @testset "dispatch disambiguator" begin
        # 3-arg ambiguous shape
        @test_throws ArgumentError find_deeper_minimum(m_sb, [1.0], [0.1])
        # 5-arg invalid shape (Minuit + x0 + errors + refit + data) mixes API styles
        @test_throws ArgumentError find_deeper_minimum(m_sb, [1.0], [0.1], refit_sb, pts_sb)
    end

    # ── refit returning wrong-length vector is treated as invalid ─────────────
    @testset "refit returning wrong-length vector is filtered" begin
        refit_bad = (subdata, start) -> fill(NaN, length(start) + 1)   # wrong length
        # All rows filtered → "only 0 valid" warning → break without adoption.
        m_out = @test_logs((:warn, r"valid resample"),
                           min_level = Logging.Warn,
                           find_deeper_minimum(m_sb, refit_bad, pts_sb;
                                               n_discovery = 4, seed = 4))
        @test m_out isa Minuit
        @test m_out.fval ≈ m_sb.fval atol = 1e-6
    end
end

@testset "find_deeper_minimum — resampling ADOPTION path (deeper basin found)" begin
    # The single-basin fixtures above only exercise the suitability-check early
    # return.  This drives the FULL adoption path: discovery → find_solution_modes
    # refine → new_min → rebuild Minuit → migrad!+hesse → loop.  It covers the
    # adoption-rebuild lines (refined_errors / prec / verify_threading / ndata).
    #
    # Tilted double well: shallow well at x[1]≈+1 (f≈+0.39), DEEP at x[1]≈−1
    # (f≈−0.41).  Start in the shallow well; `refit` returns DEEP-well candidates
    # so the refine step finds a genuinely deeper minimum and adoption fires.
    f(x) = (x[1]^2 - 1)^2 + 0.4 * x[1] + x[2]^2

    m = Minuit(f, [1.0, 0.3]; errors = [0.2, 0.2], strategy = 2)
    migrad!(m); hesse(m)
    m.ndata = 17                          # to verify ndata survives the rebuild
    @test m.values[1] > 0                 # confirm we start in the SHALLOW basin

    data = collect(1.0:10.0)
    # refit ignores the (fixed) objective's data-dependence and returns deep-well
    # candidates with a tiny bootstrap-derived jitter so the cluster is non-degenerate.
    refit_deep = (subdata, start) -> begin
        j = 0.01 * (sum(subdata) / length(subdata) - 5.5)
        return [-1.0 + j, 0.0 + j]
    end

    m_deep = find_deeper_minimum(m, refit_deep, data; n_discovery = 12, seed = 1)

    @test m_deep isa Minuit
    @test m_deep.valid                    # Minuit uses .valid, NOT is_valid
    @test m_deep.values[1] < 0            # escaped to the DEEP well
    @test m_deep.fval < m.fval - 0.5      # genuinely deeper (≈−0.41 vs ≈+0.39)
    @test m_deep.ndata == 17              # ndata carried through the adoption rebuild
    @test m_deep !== m                    # a NEW Minuit was adopted (not the input)

    # determinism: same seed ⇒ same adopted minimum
    m_deep2 = find_deeper_minimum(m, refit_deep, data; n_discovery = 12, seed = 1)
    @test collect(m_deep2.values) ≈ collect(m_deep.values)
end

@testset "find_deeper_minimum — AD-gradient parity" begin
    # Both the perturbation m::Minuit overload and the fresh-start resampling
    # overload must keep an analytical/AD gradient (not silently fall back to
    # central differences). These paths are green-but-were-untested.
    fq(x) = (x[1] - 2.0)^2 + (x[2] + 1.0)^2
    gq(x) = [2 * (x[1] - 2.0), 2 * (x[2] + 1.0)]          # exact gradient
    m_ad = Minuit(fq, [0.0, 0.0]; errors = [0.3, 0.3], grad = gq, strategy = 1)
    migrad!(m_ad)
    @test m_ad.cfwg !== nothing                           # fixture really has a gradient cf

    # ── Fix #5: perturbation m::Minuit routes through m.cfwg (AD migrad path).
    # If migrad(::CostFunctionWithGradient) / reset_ncalls! were missing this errors.
    fm = find_deeper_minimum(m_ad; n_restarts = 10, perturb = 0.4, seed = 1)
    @test JuMinuit.is_valid(fm)
    @test values(fm) ≈ [2.0, -1.0] atol = 1e-3            # converged via the AD-routed search

    # ── New fix: fresh-start resampling overload forwards cf.g for a
    # CostFunctionWithGradient.  Single-basin refit → suitability fires → returns
    # the internal m0; m0 must carry the gradient (cfwg !== nothing).
    cf_ad = m_ad.cfwg
    refit_ad = (subdata, start) -> begin
        fm2 = migrad(CostFunction(fq, 1.0), start, [0.3, 0.3]; strategy = JuMinuit.Strategy(1))
        JuMinuit.is_valid(fm2) ? collect(Float64, values(fm2)) : fill(NaN, length(start))
    end
    data_ad = fill(1.0, 16)
    # NullLogger: we only assert gradient preservation here, not the warning set
    # (an AD fit may also emit a one-time CheckGradient line at the seed).
    m_out = with_logger(NullLogger()) do
        find_deeper_minimum(cf_ad, [0.0, 0.0], [0.3, 0.3], refit_ad, data_ad;
                            n_discovery = 8, seed = 1)
    end
    @test m_out isa Minuit
    @test m_out.cfwg !== nothing                          # AD gradient preserved (the fix)

    # ── check_gradient flag preserved through the fresh-start rebuild (round-3 fix):
    # an explicit check_gradient=false must NOT be silently reset to the
    # constructor default (true).
    m_cg = Minuit(fq, [0.0, 0.0]; errors = [0.3, 0.3], grad = gq, check_gradient = false)
    @test m_cg.cfwg.check_gradient == false               # fixture really has it off
    m_out_cg = with_logger(NullLogger()) do
        find_deeper_minimum(m_cg.cfwg, [0.0, 0.0], [0.3, 0.3], refit_ad, data_ad;
                            n_discovery = 8, seed = 1)
    end
    @test m_out_cg.cfwg !== nothing
    @test m_out_cg.cfwg.check_gradient == false           # preserved, not reset to default

    # ── SITE 1 (adoption rebuild): grad + check_gradient must survive an ACTUAL
    # basin adoption, not just the suitability early-return.  Tilted double well
    # with an exact gradient; start shallow, refit returns deep-well candidates.
    fw(x) = (x[1]^2 - 1)^2 + 0.4 * x[1] + x[2]^2
    gw(x) = [4 * x[1] * (x[1]^2 - 1) + 0.4, 2 * x[2]]
    m_w = Minuit(fw, [1.0, 0.3]; errors = [0.2, 0.2], grad = gw,
                 check_gradient = false, strategy = 2)
    migrad!(m_w); hesse(m_w)
    @test m_w.values[1] > 0                               # starts in the shallow basin
    data_w = collect(1.0:10.0)
    refit_w = (subdata, start) -> begin
        j = 0.01 * (sum(subdata) / length(subdata) - 5.5)
        return [-1.0 + j, 0.0 + j]
    end
    m_w_deep = with_logger(NullLogger()) do
        find_deeper_minimum(m_w, refit_w, data_w; n_discovery = 12, seed = 1)
    end
    @test m_w_deep.values[1] < 0                          # actually adopted the deep basin
    @test m_w_deep.cfwg !== nothing                       # gradient survived adoption
    @test m_w_deep.cfwg.check_gradient == false           # check_gradient survived adoption
end
