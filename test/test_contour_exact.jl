# SPDX-License-Identifier: LGPL-2.1-or-later

@testset "contour_exact — multi-param function_cross (Phase 1.x)" begin

    @testset "function_cross_multi basic" begin
        # f(x, y, z) = (x-1)² + (y-2)² + (z-3)². Min at (1, 2, 3).
        # Fix (x, y) and scan along (1, 0) direction: minimum at z=3
        # requires no movement in z; alpha goes to root of "fval + 1".
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2 + (x[3] - 3.0)^2)
        fmin = migrad(cf, [0.0, 0.0, 0.0], [0.1, 0.1, 0.1])
        @test fmin.is_valid

        # Fix (x, y) at (1, 2) — the minimum position; ray along (1, 0).
        # The constrained 1D minimum varies as we move (x, y) away from
        # the minimum. At alpha = 1, x = 2 → fval = 1 + (y-2)² + 0 = 1.
        # So crossing at alpha = 1 (where fval = 1 = up).
        cross = JuMinuit.function_cross_multi(
            fmin, cf, [1, 2], [1.0, 2.0], [1.0, 0.0]; tlr = 0.1)
        @test cross.valid
        @test cross.aopt ≈ 1.0 atol = 0.1
    end

    @testset "contour_exact on symmetric quadratic — circle" begin
        # f(x, y) = (x-1)² + (y-2)². Minimum (1, 2), Hessian = 2·I.
        # Up = 1 → 1σ contour at radius 1 around (1, 2).
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        @test fmin.is_valid
        c = contour_exact(fmin, cf, 1, 2; npoints = 8)
        @test c.valid
        @test length(c.points) >= 4  # at least the 4 axis points
        # Every point should be at radius ≈ 1
        for (x, y) in c.points
            r = sqrt((x - 1.0)^2 + (y - 2.0)^2)
            @test r ≈ 1.0 atol = 0.15
        end
    end

    @testset "contour_exact handles correlated FCN" begin
        # f(x, y) = (x-1)² + (y-1)² + 0.5·x·y. Hessian has off-diagonal.
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 1.0)^2 + 0.5 * x[1] * x[2])
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        @test fmin.is_valid
        c = contour_exact(fmin, cf, 1, 2; npoints = 8)
        @test c.valid
        # Boundary points should all be at reasonable distance from min
        center = Base.values(fmin)
        for (x, y) in c.points
            d = sqrt((x - center[1])^2 + (y - center[2])^2)
            @test 0.2 < d < 10.0
        end
    end

    @testset "Argument validation" begin
        cf = CostFunction(x -> sum(abs2, x))
        fmin = migrad(cf, [1.0, 2.0], [0.1, 0.1])
        @test_throws ArgumentError contour_exact(fmin, cf, 0, 2)
        @test_throws ArgumentError contour_exact(fmin, cf, 1, 3)
        @test_throws ArgumentError contour_exact(fmin, cf, 1, 1)
        @test_throws ArgumentError contour_exact(fmin, cf, 1, 2; npoints = 3)
    end

    @testset "function_cross_multi argument validation" begin
        cf = CostFunction(x -> sum(abs2, x))
        fmin = migrad(cf, [1.0, 2.0], [0.1, 0.1])
        @test_throws DimensionMismatch JuMinuit.function_cross_multi(
            fmin, cf, [1], [1.0, 2.0], [1.0])
        @test_throws DimensionMismatch JuMinuit.function_cross_multi(
            fmin, cf, [1, 2], [1.0, 2.0], [1.0])
        # n == npar (no free parameters) is now supported via the
        # all-fixed degenerate path used by 2D contour.
        cr = JuMinuit.function_cross_multi(
            fmin, cf, [1, 2], [1.0, 2.0], [1.0, 0.0])
        @test cr.aopt isa Float64 || isnan(cr.aopt)
    end
end
