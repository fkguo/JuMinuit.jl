# SPDX-License-Identifier: LGPL-2.1-or-later

using LinearAlgebra: Symmetric

@testset "covariance_squeeze.jl — MnCovarianceSqueeze" begin
    prec = MachinePrecision()

    @testset "squeeze_symmetric — pure drop row/col" begin
        # Construct a 3×3 symmetric matrix and drop the middle row/col
        M = Float64[1 2 3; 0 4 5; 0 0 6]
        S = Symmetric(M, :U)
        Sq = JuMinuit.squeeze_symmetric(S, 2)
        @test size(Sq) == (2, 2)
        # Remaining entries: (1,1), (1,3), (3,1)=(1,3), (3,3) → new indices (1,1), (1,2), (2,1), (2,2)
        @test Sq[1, 1] == 1.0
        @test Sq[1, 2] == 3.0
        @test Sq[2, 2] == 6.0
        @test Sq[2, 1] == 3.0  # symmetric read

        # Drop first row/col
        Sq1 = JuMinuit.squeeze_symmetric(S, 1)
        @test Sq1[1, 1] == 4.0
        @test Sq1[1, 2] == 5.0
        @test Sq1[2, 2] == 6.0

        # Drop last row/col
        Sq3 = JuMinuit.squeeze_symmetric(S, 3)
        @test Sq3[1, 1] == 1.0
        @test Sq3[1, 2] == 2.0
        @test Sq3[2, 2] == 4.0
    end

    @testset "squeeze_symmetric — error cases" begin
        S = Symmetric([1.0 2.0; 0.0 3.0], :U)
        @test_throws ArgumentError JuMinuit.squeeze_symmetric(S, 0)
        @test_throws ArgumentError JuMinuit.squeeze_symmetric(S, 3)
        # 1x1 cannot squeeze
        S1 = Symmetric(reshape([1.0], 1, 1), :U)
        @test_throws ArgumentError JuMinuit.squeeze_symmetric(S1, 1)
    end

    @testset "squeeze_error — invert/squeeze/re-invert round trip" begin
        # Set up a 3×3 well-conditioned inverse-Hessian (V).
        # H = inv(V) is the Hessian. Squeezing param 2 from V should
        # produce a 2×2 V' equal to inv(squeeze(H, 2)).
        V = Symmetric(Float64[
            2.0  0.1  0.05;
            0.1  3.0  0.2;
            0.05 0.2  1.5
        ], :U)
        err = MinimumError(V, 0.001)
        sq = JuMinuit.squeeze_error(err, 2; prec)
        @test size(sq) == (2, 2)
        @test is_valid(sq)
        @test sq.dcovar == 0.001  # preserved per C++

        # Independent oracle: build H = inv(V), drop row/col 2, invert back
        H_full = inv(Matrix(V))
        H_squeezed = H_full[[1, 3], [1, 3]]
        V_expected = inv(H_squeezed)
        for i in 1:2, j in 1:2
            @test sq.inv_hessian[i, j] ≈ V_expected[i, j] atol = 1e-10
        end
    end

    @testset "squeeze_error — singular matrix triggers diagonal fallback" begin
        # Construct a rank-1 inverse-Hessian (will fail to invert)
        v = [1.0, 2.0, 3.0]
        # V = v·v' is rank-1 (singular); sym_invert! will throw
        Vsing_mat = v * v'
        # Make sure the construction is positive-semidefinite Symmetric
        Vsing = Symmetric(0.5 * (Vsing_mat + Vsing_mat'), :U)
        err = MinimumError(Vsing, 0.5)
        sq = JuMinuit.squeeze_error(err, 2; prec)
        @test invert_failed(sq)
        @test size(sq) == (2, 2)
    end
end
