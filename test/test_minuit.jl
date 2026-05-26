# SPDX-License-Identifier: LGPL-2.1-or-later

@testset "minuit.jl — iminuit-style Minuit wrapper" begin

    @testset "Constructor + property access (no migrad)" begin
        m = Minuit(x -> sum(abs2, x), [1.0, 2.0];
                    names = ["x", "y"], errors = [0.1, 0.2])
        @test m isa Minuit
        @test m.ndim == 2
        @test m.npar == 2
        @test m.values == [1.0, 2.0]  # initial
        @test m.errors == [0.1, 0.2]  # initial
        @test isnan(m.fval)
        @test isnan(m.edm)
        @test m.nfcn == 0
        @test !m.valid
        @test m.covariance === nothing
    end

    @testset "migrad! workflow" begin
        m = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2, [0.0, 0.0];
                    names = ["x", "y"], errors = [0.1, 0.1])
        migrad!(m)
        @test m.valid
        @test m.fval < 1e-8
        @test m.values[1] ≈ 1.0 atol = 1e-4
        @test m.values[2] ≈ 2.0 atol = 1e-4
        @test m.covariance isa Matrix{Float64}
        @test size(m.covariance) == (2, 2)
    end

    @testset "Bounded + fixed parameters via Minuit" begin
        m = Minuit(x -> (x[1] - 0.5)^2 + (x[2] - 3.0)^2, [0.3, 5.0];
                    names = ["a", "b"], errors = [0.1, 0.1],
                    limits = [(0.0, 1.0), nothing],
                    fixed = [false, true])
        migrad!(m)
        @test m.valid
        @test m.values[1] ≈ 0.5 atol = 0.01
        @test m.values[2] == 5.0  # fixed bit-exact
        @test m.errors[2] == 0.0   # fixed → no error
    end

    @testset "minos! workflow" begin
        m = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2, [0.0, 0.0];
                    names = ["x", "y"])
        migrad!(m)
        minos!(m, 1)
        @test haskey(m.minos_errors, 1)
        e = m.minos_errors[1]
        @test JuMinuit.is_valid(e)
        @test e.upper ≈ 1.0 atol = 0.1
        @test e.lower ≈ -1.0 atol = 0.1

        # By name
        minos!(m, "y")
        @test haskey(m.minos_errors, 2)

        # All free
        m2 = Minuit(x -> sum(abs2, x), [1.0, 2.0])
        migrad!(m2)
        minos!(m2)
        @test length(m2.minos_errors) == 2
    end

    @testset "contour workflow" begin
        m = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2, [0.0, 0.0];
                    names = ["x", "y"])
        migrad!(m)
        c = contour(m, 1, 2; npoints = 10)
        @test c isa ContoursError
        @test c.valid
        @test length(c.points) == 10

        # By name
        c2 = contour(m, "x", "y"; npoints = 8)
        @test c2.valid
        @test length(c2.points) == 8
    end

    @testset "Pretty print" begin
        m = Minuit(x -> sum(abs2, x), [1.0, 2.0];
                    names = ["x", "y"],
                    limits = [(-5.0, 5.0), nothing])
        # Before migrad
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), m)
        s = String(take!(buf))
        @test occursin("not yet minimized", s)
        @test occursin("[-5.0, 5.0]", s)

        # After migrad — Phase 3 C1 Unicode table format
        migrad!(m)
        buf2 = IOBuffer()
        show(buf2, MIME"text/plain"(), m)
        s2 = String(take!(buf2))
        # Header line carries fval / edm / nfcn / status
        @test occursin("fval=", s2)
        @test occursin("nfcn=", s2)
        @test occursin("Valid", s2)
        # Unicode box-drawing characters present
        @test occursin("┌", s2)
        @test occursin("┤", s2)
        @test occursin("└", s2)
        # Column headers
        for col in ("Name", "Value", "Hesse ±", "Minos −", "Minos +",
                    "Limit −", "Limit +", "Fixed")
            @test occursin(col, s2)
        end
    end

    @testset "C1 (a) at-limit warning detection" begin
        # Force a parameter to sit on a tight lower bound: fit
        # (x-0.5)² with x ∈ [0.3, 10]. The minimum is at 0.5, the lower
        # bound is 0.3 away — well within 1σ (Hesse err ≈ 1.0).
        cf = x -> (x[1] - 0.5)^2
        m = Minuit(cf, [1.0]; name = ["a"], limit_a = (0.3, 10.0))
        migrad(m)
        @test m.is_valid
        # At-limit detector should flag `a` (lower edge within 1σ)
        @test 1 in JuMinuit._at_limit_indices(m)
        # Warning visible in text/plain output
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), m)
        s = String(take!(buf))
        @test occursin("⚠", s)
        @test occursin("`a`", s)
        @test occursin("lower limit", s)
        @test occursin("unreliable", s)

        # Negative case: parameter NOT at limit
        cf2 = x -> (x[1] - 5.0)^2
        m2 = Minuit(cf2, [4.0]; name = ["a"], limit_a = (-100.0, 100.0))
        migrad(m2)
        @test isempty(JuMinuit._at_limit_indices(m2))
        buf2 = IOBuffer()
        show(buf2, MIME"text/plain"(), m2)
        s2 = String(take!(buf2))
        @test !occursin("⚠", s2)
    end

    @testset "C1 (c) HTML repr (IJulia / Pluto)" begin
        cf = x -> sum(abs2, x .- [1.0, 2.0])
        m = Minuit(cf, [0.0, 0.0]; names = ["a", "b"])
        # Before migrad
        buf = IOBuffer()
        show(buf, MIME"text/html"(), m)
        @test occursin("not yet minimized", String(take!(buf)))

        # After migrad
        migrad(m)
        buf = IOBuffer()
        show(buf, MIME"text/html"(), m)
        s = String(take!(buf))
        # HTML structure
        @test occursin("<table", s)
        @test occursin("<thead", s)
        @test occursin("<tbody", s)
        @test occursin("</table>", s)
        # Column headers and a status badge
        @test occursin("Hesse ±", s)
        @test occursin("Minos −", s)
        @test occursin("Valid", s)

        # at-limit + HTML: yellow warning div appears
        m_lim = Minuit(x -> (x[1] - 0.5)^2, [1.0];
                       name = ["a"], limit_a = (0.3, 10.0))
        migrad(m_lim)
        buf = IOBuffer()
        show(buf, MIME"text/html"(), m_lim)
        s_lim = String(take!(buf))
        @test occursin("⚠", s_lim)
        @test occursin("<code>a</code>", s_lim)
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError Minuit(x -> 0.0, [1.0, 2.0]; names = ["x"])
        m = Minuit(x -> sum(abs2, x), [1.0, 2.0])
        @test_throws ArgumentError minos!(m, 1)  # no migrad! yet
        @test_throws ArgumentError JuMinuit.hesse(m)  # no migrad! yet
    end

    @testset "hesse(m) refreshes the covariance (was placeholder)" begin
        # Task #36 — `hesse(m::Minuit)` used to be a no-op. After the
        # fix, a Strategy(0) MIGRAD followed by `hesse(m)` should give
        # a covariance numerically close to a Strategy(2) MIGRAD's
        # output (both end with a full numerical-HESSE pass).
        cf_fn = x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2
        m_s0 = Minuit(cf_fn, [0.0, 0.0]; errors = [0.1, 0.1])
        migrad(m_s0; strategy = Strategy(0))
        cov_s0_dfp = collect(m_s0.covariance)   # DFP estimate

        JuMinuit.hesse(m_s0; strategy = Strategy(1))
        cov_s0_hesse = collect(m_s0.covariance)   # numerical HESSE
        @test m_s0.is_valid
        # Strategy(2) MIGRAD also ends with numerical HESSE.
        m_s2 = Minuit(cf_fn, [0.0, 0.0]; errors = [0.1, 0.1])
        migrad(m_s2; strategy = Strategy(2))
        cov_s2 = collect(m_s2.covariance)
        @test m_s2.is_valid

        # The hesse(m) refresh should match Strategy(2) MIGRAD's cov
        # element-by-element (both compute numerical 2nd-derivative
        # Hessian at the same converged minimum). Pure quadratic FCN
        # so the inverse Hessian is exact in both paths.
        @test cov_s0_hesse ≈ cov_s2 atol = 1e-8

        # Bounded variant: ensure the int↔ext transform round-trips
        # in hesse(m) too.
        m_bnd = Minuit(cf_fn, [0.0, 0.0]; errors = [0.1, 0.1],
                                          limit_x0 = (-5.0, 5.0))
        migrad(m_bnd; strategy = Strategy(0))
        cov_bnd_pre = collect(m_bnd.covariance)
        JuMinuit.hesse(m_bnd; strategy = Strategy(1))
        @test m_bnd.is_valid
        cov_bnd_post = collect(m_bnd.covariance)
        # Bounded path: covariance shape preserved, diagonals positive.
        @test size(cov_bnd_post) == size(cov_bnd_pre)
        @test all(diag(cov_bnd_post) .> 0)
        # External errors get refreshed via int2ext_error.
        @test all(m_bnd.errors .> 0)

        # hesse(m) returns m for chaining.
        m_chain = Minuit(cf_fn, [0.0, 0.0]; errors = [0.1, 0.1])
        migrad(m_chain)
        @test JuMinuit.hesse(m_chain) === m_chain
    end
end
