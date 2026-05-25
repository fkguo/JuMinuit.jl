# SPDX-License-Identifier: LGPL-2.1-or-later

@testset "hesse.jl — MnHesse" begin

    @testset "Quadratic — exact Hessian recovery" begin
        # f(x) = Σ xᵢ². Hessian = 2·I. Inverse = 0.5·I.
        cf = CostFunction(x -> sum(abs2, x))
        # Run MIGRAD to convergence first
        m = migrad(cf, [1.0, 2.0, 3.0], [0.1, 0.1, 0.1])
        @test m.is_valid
        # Apply HESSE
        new_state = hesse(cf, m.state, Strategy(1))
        @test is_valid(new_state.error)
        @test is_accurate(new_state.error)  # dcov < 0.1
        # The inv_hessian should be ≈ 0.5·I for the quadratic
        for i in 1:3, j in 1:3
            expected = i == j ? 0.5 : 0.0
            @test new_state.error.inv_hessian[i, j] ≈ expected atol = 1e-5
        end
    end

    @testset "Hesse refines off-diagonal" begin
        # f(x, y) = x² + y² + 0.1·x·y; Hessian = [2 0.1; 0.1 2].
        # Inverse = (1/(4-0.01)) · [2 -0.1; -0.1 2] ≈ 0.5006·[2 -0.1; -0.1 2]
        cf = CostFunction(x -> x[1]^2 + x[2]^2 + 0.1 * x[1] * x[2])
        m = migrad(cf, [1.0, 1.0], [0.1, 0.1])
        @test m.is_valid
        st = hesse(cf, m.state, Strategy(1))
        # 2x2 Hessian: H = [2 0.1; 0.1 2]; det = 3.99; inv = (1/3.99)·[2 -0.1; -0.1 2]
        inv_det = 1.0 / (4.0 - 0.01)
        @test st.error.inv_hessian[1, 1] ≈ inv_det * 2.0 atol = 1e-4
        @test st.error.inv_hessian[2, 2] ≈ inv_det * 2.0 atol = 1e-4
        @test st.error.inv_hessian[1, 2] ≈ -inv_det * 0.1 atol = 1e-4
    end

    @testset "Hesse on 1D quadratic" begin
        cf = CostFunction(x -> 3.0 * x[1]^2)  # Hessian = 6; inv = 1/6
        m = migrad(cf, [1.0], [0.1])
        st = hesse(cf, m.state, Strategy(1))
        @test st.error.inv_hessian[1, 1] ≈ 1 / 6 atol = 1e-5
        @test is_valid(st.error)
    end

    @testset "Hesse increments NFcn" begin
        cf = CostFunction(x -> sum(abs2, x))
        m = migrad(cf, [1.0, 2.0], [0.1, 0.1])
        nfcn_before = nfcn(m.state)
        st = hesse(cf, m.state, Strategy(1))
        # HESSE does at least 2n + n·(n-1)/2 + 1 calls = 2·2 + 1 + 1 = 6
        @test nfcn(st) >= nfcn_before + 2
    end

    @testset "Strategy levels exercise different ncycles" begin
        cf = CostFunction(x -> sum(abs2, x))
        m = migrad(cf, [1.0, 1.0], [0.1, 0.1])
        # Strategy 0/1/2 differ in hessian_ncycles (3/5/7); all should
        # converge on a smooth quadratic.
        for level in (0, 1, 2)
            st = hesse(cf, m.state, Strategy(level))
            @test is_valid(st.error)
        end
    end

    @testset "Type stability" begin
        cf = CostFunction(x -> sum(abs2, x))
        m = migrad(cf, [1.0, 2.0], [0.1, 0.1])
        @test (@inferred hesse(cf, m.state, Strategy(1))) isa MinimumState
    end
end
