# SPDX-License-Identifier: LGPL-2.1-or-later

@testset "plot_text.jl — mn_plot_text ASCII renderer (gap M2)" begin

    # Unicode box-drawing characters used by the frame. Tests
    # rendering of these reuses src/minuit.jl's pretty-print palette.
    box_chars = ('┌', '─', '┐', '│', '└', '┘')

    @testset "circle (20 points): basic invariants" begin
        npts = 20
        θs = range(0, 2π; length = npts + 1)[1:npts]
        pts = [(cos(θ), sin(θ)) for θ in θs]

        out = mn_plot_text(pts; width = 50, height = 18,
                            par_x = "alpha", par_y = "beta")

        @test out isa String
        @test !isempty(out)
        # > 1 line — header (3 lines) + top frame + body + bottom + legend ≥ 7 lines
        nlines = count(==('\n'), out)
        @test nlines > 5

        # Parameter names appear in the header
        @test occursin("alpha", out)
        @test occursin("beta", out)
        @test occursin("Δx", out)  # scale legend on header line

        # Box-drawing frame present
        for c in box_chars
            @test occursin(c, out)
        end

        # All 20 sample points reported
        @test occursin("(20 points)", out)
        # No "EMPTY" sentinel
        @test !occursin("EMPTY", out)
    end

    @testset "ContoursError dispatch" begin
        cf = CostFunction(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2)
        fmin = migrad(cf, [0.0, 0.0], [0.1, 0.1])
        c = contour(fmin, cf, 1, 2; npoints = 16)
        @test c.valid

        out = mn_plot_text(c; width = 50, height = 18)
        @test out isa String
        # Recipe-style "par[1]" / "par[2]" header labels
        @test occursin("par[1]", out)
        @test occursin("par[2]", out)
        # Minimum should be marked with X (centroid drawn from MINOS)
        @test occursin('X', out)
        @test occursin("X = minimum", out)
    end

    @testset "empty pts → graceful EMPTY message" begin
        out = mn_plot_text(Tuple{Float64,Float64}[]; width = 40, height = 10)
        @test out isa String
        @test occursin("EMPTY", out)
        # Must not crash, must not contain frame characters
        @test !occursin('┌', out)
    end

    @testset "all-non-finite → behaves like empty" begin
        pts = [(NaN, 1.0), (2.0, Inf), (-Inf, 0.0)]
        out = mn_plot_text(pts)
        @test occursin("EMPTY", out)
    end

    @testset "mixed finite + non-finite: non-finite skipped" begin
        pts = [(1.0, 1.0), (NaN, 2.0), (2.0, Inf), (3.0, 3.0)]
        out = mn_plot_text(pts; width = 30, height = 10)
        # Only 2 finite points — rest silently dropped
        @test occursin("(2 points)", out)
        @test !occursin("EMPTY", out)
    end

    @testset "width / height kwargs change line count" begin
        npts = 12
        θs = range(0, 2π; length = npts + 1)[1:npts]
        pts = [(cos(θ), sin(θ)) for θ in θs]

        small = mn_plot_text(pts; width = 20, height = 10)
        large = mn_plot_text(pts; width = 80, height = 40)
        n_small = count(==('\n'), small)
        n_large = count(==('\n'), large)
        @test n_large > n_small
        # And the wide form should also produce visibly wider top-frame line
        top_small = split(small, '\n')[4]   # 3-line header + frame line
        top_large = split(large, '\n')[4]
        @test length(top_large) > length(top_small)
    end

    @testset "overprint: differing characters → `&`" begin
        # Per C++ mnplot semantics (mntplot.cxx:162) two SAME-character
        # stamps coalesce silently — only differing characters collide
        # to '&'. Verify both halves of the contract.
        same_char = mn_plot_text([(0.5, 0.5), (0.5, 0.5)];
                                   width = 20, height = 10)
        @test !occursin('&', same_char)

        # '*' at (0.5, 0.5) + 'X' centroid at same cell → '&' fires.
        out = mn_plot_text([(0.5, 0.5)]; width = 20, height = 10,
                             x_center = (0.5, 0.5))
        @test occursin('&', out)
        @test occursin("overlap", out)
    end

    @testset "degenerate single-point input does not crash" begin
        out = mn_plot_text([(1.0, 1.0)]; width = 20, height = 10)
        @test out isa String
        @test occursin("(1 points)", out)
        # Frame must still be present
        @test occursin('┌', out)
        @test occursin('└', out)
    end

    @testset "_mn_bins smoke checks (round numbers)" begin
        # Round binning per Minuit2 mnbins: width-10 → expects ≈ 1.0/0.5/0.2/0.1 spacing
        bl, bh, nb, bwid = JuMinuit._mn_bins(-1.0, 1.0, 10)
        @test bl <= -1.0 && bh >= 1.0
        @test nb >= 1
        @test bwid > 0
        # The width times bin count must enclose the range
        @test bl + nb * bwid ≈ bh atol=1e-9
    end
end
