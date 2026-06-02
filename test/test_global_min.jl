# SPDX-License-Identifier: LGPL-2.1-or-later
using JuMinuit
using Test

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
