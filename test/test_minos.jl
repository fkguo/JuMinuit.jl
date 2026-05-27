# SPDX-License-Identifier: LGPL-2.1-or-later

@testset "minos.jl + function_cross.jl" begin

    @testset "MinosError struct" begin
        # min_par_value field stores the parameter value at the minimum,
        # mirroring C++ MinosError::Min() (parallel-review #4 B2 fix).
        e = MinosError(1, 1.5, 0.5, -0.5, true, true, false, false, false, false, 100)
        @test e.par_idx == 1
        @test e.min_par_value == 1.5
        @test e.upper == 0.5
        @test e.lower == -0.5
        @test JuMinuit.is_valid(e)
    end

    @testset "Symmetric quadratic — MINOS ≈ Hesse" begin
        # f(x, y) = (x - 1)² + (y - 2)². Minimum at (1, 2), fval = 0.
        # Hessian is 2·I, so V = 0.5·I, errors = sqrt(2·1·0.5) = 1.0.
        # MINOS should give upper = -lower ≈ 1.0 for each parameter.
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        @test fmin.is_valid

        e1 = minos(fmin, cf, 1)
        @test JuMinuit.is_valid(e1)
        # Symmetric → upper ≈ -lower
        @test e1.upper ≈ 1.0 atol = 0.1
        @test e1.lower ≈ -1.0 atol = 0.1

        e2 = minos(fmin, cf, 2)
        @test JuMinuit.is_valid(e2)
        @test e2.upper ≈ 1.0 atol = 0.1
        @test e2.lower ≈ -1.0 atol = 0.1
    end

    @testset "All-parameters convenience" begin
        cf = CostFunction(x -> sum(abs2, x .- [1.0, 2.0]))
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        errs = minos(fmin, cf)
        @test length(errs) == 2
        @test all(JuMinuit.is_valid, errs)
    end

    @testset "Argument validation" begin
        cf = CostFunction(x -> sum(abs2, x))
        fmin = migrad(cf, [1.0, 2.0], [0.1, 0.1])
        @test_throws ArgumentError minos(fmin, cf, 0)
        @test_throws ArgumentError minos(fmin, cf, 3)

        # n=1 should throw (cannot fix the only parameter)
        cf1 = CostFunction(x -> x[1]^2)
        fmin1 = migrad(cf1, [1.0], [0.1])
        @test_throws ArgumentError minos(fmin1, cf1, 1)
    end

    @testset "function_cross — direct call" begin
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        # Upper direction along param 1
        cr_up = JuMinuit.function_cross(fmin, cf, 1, +1.0)
        @test cr_up.valid
        @test cr_up.aopt > 0.5
        # Lower direction along param 1
        cr_lo = JuMinuit.function_cross(fmin, cf, 1, -1.0)
        @test cr_lo.valid
        @test cr_lo.aopt > 0.5  # aopt is the magnitude regardless of sign
    end

    @testset "function_cross — parabolic path (A3/A4)" begin
        # Phase 1.x A3/A4 (parallel-review #4) — non-quadratic CF that
        # exercises the L500 MnParabola 3-point fit. With a quartic
        # term in x[1], the crossing surface is f = a·(x-1)⁴ + (y-2)²,
        # so the level set at f = fmin+1 is x = 1 ± 1/a^(1/4). The
        # crossing α (relative to the post-fit MIGRAD step σ_x) should
        # be ~1·σ_x and the parabolic fit converges in one or two L500
        # iterations vs many for the linear-only path.
        cf = CostFunction(x -> 4.0 * (x[1] - 1.0)^4 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        @test JuMinuit.is_valid(fmin)
        cr_up = JuMinuit.function_cross(fmin, cf, 1, +1.0)
        @test cr_up.valid
        # For x[1], the crossing is at x = 1 + (1/4)^(1/4) ≈ 1.707;
        # the 1σ post-fit step is also nonquadratic-skewed but the
        # parabolic search still converges (rough magnitude check
        # only; exact α depends on σ_x from the converged Hessian).
        @test cr_up.aopt > 0.0
        @test cr_up.nfcn < 1500  # parabola fit should NOT explode call count
    end

    @testset "parabola helpers — direct unit tests" begin
        # A·x² + B·x + C through (0, 1), (1, 0), (2, 1)
        # Expected: A=1, B=-2, C=1   (i.e., f(x) = (x-1)²)
        A, B, C = JuMinuit._parabola_fit3([0.0, 1.0, 2.0], [1.0, 0.0, 1.0])
        @test A ≈ 1.0
        @test B ≈ -2.0
        @test C ≈ 1.0

        # Solve (x-1)² = 2 → roots 1 ± √2. Positive-slope root is 1+√2.
        prec = JuMinuit.MachinePrecision()
        sol = JuMinuit._parabola_solve_for_aim(1.0, -2.0, 1.0, 2.0, prec)
        @test sol !== nothing
        x_sol, slope = sol
        @test x_sol ≈ 1.0 + sqrt(2.0)
        @test slope > 0  # positive-slope root selected

        # Negative-curvature (A < 0) parabola: f(x) = -(x-1)² + 2.
        # Solve = 1 needs determ = B² - 4A(C-aim) = 4 - 4·(-1)·(2-1-2) = 4 - 4 = 0
        # → single root x=1. Discriminant ≥ 0 means we still get a result.
        sol2 = JuMinuit._parabola_solve_for_aim(-1.0, 2.0, 1.0, 1.0, prec)
        @test sol2 !== nothing
        # Negative curvature with too-high aim → discriminant < 0
        sol3 = JuMinuit._parabola_solve_for_aim(-1.0, 2.0, 1.0, 5.0, prec)
        @test sol3 === nothing
    end

    @testset "three-point classifier — direct unit tests" begin
        # 3 points around aim=0: f = (-1, -0.5, +1). noless=2, ibest=2 (closest to 0).
        ibest, iworst, ileft, iright, iout, noless, ecmn, ecmx =
            JuMinuit._three_point_classify([0.0, 0.5, 1.0], [-1.0, -0.5, 1.0], 0.0)
        @test noless == 2
        @test ibest == 2          # |−0.5−0| = 0.5 is smallest
        @test iworst == 1         # |−1−0| = 1, |1−0| = 1; first-seen wins iworst
        @test iright == 3         # f[3] = 1 > 0 → right side
        @test ileft == 2          # ileft tracks closest-to-aim on left; f[2]=-0.5 > f[1]=-1
        @test iout == 1           # the farther-left point becomes redundant

        # All three above aim: noless=0
        _, _, _, _, _, noless0, _, _ =
            JuMinuit._three_point_classify([0.0, 1.0, 2.0], [2.0, 3.0, 5.0], 1.0)
        @test noless0 == 0

        # All three below aim: noless=3
        _, _, _, _, _, noless3, _, _ =
            JuMinuit._three_point_classify([0.0, 1.0, 2.0], [-2.0, -1.0, -0.5], 1.0)
        @test noless3 == 3

        # default_ibest tie-break (Opus review IMPORTANT #5): when all three
        # |f - aim| are equal, the initial classifier uses default_ibest=3.
        ib_init, _, _, _, _, _, _, _ =
            JuMinuit._three_point_classify([0.0, 0.5, 1.0], [1.0, 1.0, 1.0], 0.0;
                                            default_ibest = 3)
        @test ib_init == 3
        # L500 classifier uses default_ibest=1
        ib_l500, _, _, _, _, _, _, _ =
            JuMinuit._three_point_classify([0.0, 0.5, 1.0], [1.0, 1.0, 1.0], 0.0;
                                            default_ibest = 1)
        @test ib_l500 == 1
    end

    @testset "function_cross — tlr=0.01 crossing tolerance (BLOCKING #1)" begin
        # Opus review BLOCKING #1 — C++ MnFunctionCross.cxx:38-40 hardcodes
        # the CROSSING convergence to tlr=0.01 regardless of user-supplied
        # tlr (which only controls inner-MIGRAD via 0.5·tlr). Earlier the
        # Julia code propagated user tlr (default 0.1) to the convergence
        # check too → 10× looser than C++. The fix should give aopt within
        # ~1% of the analytic answer even at the loose default.
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        # The analytic 1σ crossing along x[1] is at α = 1.0 (since
        # σ_x = sqrt(2·1·0.5) = 1 and the crossing at f=fmin+1 is x = 1+1·σ_x).
        cr = JuMinuit.function_cross(fmin, cf, 1, +1.0)
        @test cr.valid
        # Tight 1% tolerance — without the C++ tlr=0.01 override this
        # would have a 10× looser allowable error.
        @test cr.aopt ≈ 1.0 atol = 0.01

        # Even when the user passes a loose tlr=0.5 (deliberately huge),
        # the crossing tlf=0.01·up should still pin aopt within ~1%
        # because the override decouples user-tlr from the convergence
        # check (only inner-MIGRAD sees `tol = 0.5·tlr`).
        cr_loose = JuMinuit.function_cross(fmin, cf, 1, +1.0; tlr = 0.5)
        @test cr_loose.valid
        @test cr_loose.aopt ≈ 1.0 atol = 0.05  # inner-MIGRAD looseness only
    end

    @testset "function_cross — non-quadratic many-iterations (BLOCKING #2)" begin
        # Opus review BLOCKING #2 — the C++ fall-through branch
        # (`alsb[iworst] = alsb[2]; goto L460`) is hit when the third
        # probe lands closer to aim than the first two but all 3 stay
        # one-sided. This is common on non-quadratic level surfaces
        # where the initial linear extrapolation overshoots. Without
        # the fall-through, the algorithm would have returned invalid.
        # Verify that on a quartic CF (non-quadratic crossing), the
        # algorithm DOES converge to a valid crossing.
        cf = CostFunction(x -> 4.0 * (x[1] - 1.0)^4 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        cr = JuMinuit.function_cross(fmin, cf, 1, +1.0)
        @test cr.valid                # would be false if fall-through missing
        @test cr.aopt > 0.0
        @test cr.nfcn < 1500          # bounded call count
    end

    @testset "_fix_one_param / _fix_multi_params — zero per-call alloc (V3 lift)" begin
        # Phase A V3 — perf-regression guards. The fix-* wrappers MUST NOT
        # allocate per call (lifted full_buf in the closure). If a future
        # refactor reintroduces the per-call `Vector{Float64}(undef, n_)`
        # alloc, these tests catch it immediately. 18% wall-time win on
        # all corpus benchmarks (rosenbrock_10d / gauss_ll_10_1000 / quad_4d)
        # depends on the zero-alloc invariant.
        cf = CostFunction(x -> sum(abs2, x), 1.0)
        cf_one = JuMinuit._fix_one_param(cf, 3, 0.5, 5)
        y4 = [0.1, 0.2, 0.3, 0.4]
        # Warmup (compile)
        cf_one(y4)
        # Two consecutive calls must both be zero-alloc — guards against
        # accidental closure repromotion under future precompile changes.
        @test (@allocated cf_one(y4)) == 0
        @test (@allocated cf_one(y4)) == 0
        # Return-type stability (the wrapped FCN returns Float64 → wrapper
        # must too; @inferred fails if Julia infers Any/Union).
        @test (@inferred cf_one(y4)) isa Float64

        cf_multi = JuMinuit._fix_multi_params(cf, [1, 3], [0.5, 0.5], 5)
        y3 = [0.1, 0.2, 0.3]
        cf_multi(y3)
        @test (@allocated cf_multi(y3)) == 0
        @test (@allocated cf_multi(y3)) == 0
        @test (@inferred cf_multi(y3)) isa Float64

        # Numerical-correctness sanity: splicing fixed + free params produces
        # the same value as a manual splice — guards against off-by-one in
        # the lifted-buffer write pattern.
        # cf_one: par at index 3 fixed to 0.5; free = [0.1, 0.2, 0.3, 0.4]
        @test cf_one(y4) ≈ 0.1^2 + 0.2^2 + 0.5^2 + 0.3^2 + 0.4^2
        # cf_multi: par at indices 1, 3 fixed to 0.5, 0.5; free = [0.1,0.2,0.3]
        @test cf_multi(y3) ≈ 0.5^2 + 0.1^2 + 0.5^2 + 0.2^2 + 0.3^2
    end
end
