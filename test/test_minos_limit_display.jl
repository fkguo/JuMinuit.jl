# SPDX-License-Identifier: LGPL-2.1-or-later

using JSON
using NativeMinuit
using RecipesBase
using Test

@testset "MINOS at-limit public presentation" begin
    min_value = 3.7799991402523414
    upper_distance = 3.78 - min_value

    function merror(; upper = 0.2, lower = -0.3,
                    upper_valid = true, lower_valid = true,
                    upper_par_limit = false, lower_par_limit = false)
        return MinosError(
            1, min_value, upper, lower,
            upper_valid, lower_valid,
            false, false, false, false,
            upper_par_limit, lower_par_limit,
            17,
        )
    end

    upper_limit = merror(upper = upper_distance, upper_par_limit = true)
    lower_limit = merror(lower = -0.15, lower_par_limit = true)
    both_limits = merror(
        upper = upper_distance,
        lower = -0.15,
        upper_par_limit = true,
        lower_par_limit = true,
    )
    crossing = merror()
    upper_invalid = merror(upper = 0.4, upper_valid = false)
    zero_upper = merror(upper = 0.0)
    both_invalid = merror(upper_valid = false, lower_valid = false)
    upper_limit_lower_invalid = merror(
        upper = upper_distance,
        upper_par_limit = true,
        lower_valid = false,
    )
    large_upper_limit = merror(
        upper = 1.23456789e100,
        upper_par_limit = true,
    )

    function bounded_profile(; curvature, center = 0.0, start = 0.2,
                             limits = (-10.0, 10.0))
        f = p -> curvature * (p[1] - center)^2 + (p[2] - 0.4)^2
        m = Minuit(
            f, [start, 0.0];
            names = ["x", "nuisance"],
            errors = [0.1, 0.1],
            limits = [limits, nothing],
        )
        migrad!(m)
        minos!(m, 1)
        return m, m.minos_errors[1]
    end

    @testset "real bounded profile regressions" begin
        @testset "wide quadratic finds both confidence crossings" begin
            m, e = bounded_profile(curvature = 1.0, limits = (-10.0, 10.0))
            @test e.upper_valid && e.lower_valid
            @test !e.upper_par_limit && !e.lower_par_limit
            @test NativeMinuit.has_closed_interval(e)

            single_plain = sprint(show, MIME"text/plain"(), e)
            table_plain = sprint(show, MIME"text/plain"(), m)
            table_html = sprint(show, MIME"text/html"(), m)
            table_tex = to_latex(m)
            @test occursin("MINOS error", single_plain)
            for rendered in (single_plain, table_plain, table_html, table_tex)
                @test !occursin("boundary displacement", lowercase(rendered))
                @test !occursin("scan reached limit", lowercase(rendered))
            end
            @test occursin("+", table_plain) && occursin("−", table_plain)
            @test occursin("<sup>+", table_html) && occursin("<sub>−", table_html)
            @test occursin("^{+", table_tex) && occursin("_{-", table_tex)
        end

        @testset "interior minimum with flat profile reaches both limits" begin
            m, e = bounded_profile(curvature = 0.01, limits = (-1.0, 1.0))
            @test -1.0 < m.values[1] < 1.0
            @test e.upper_valid && e.lower_valid
            @test e.upper_par_limit && e.lower_par_limit
            @test is_valid(e)
            @test !NativeMinuit.has_closed_interval(e)

            single_plain = lowercase(sprint(show, MIME"text/plain"(), e))
            table_plain = lowercase(sprint(show, MIME"text/plain"(), m))
            table_html = lowercase(sprint(show, MIME"text/html"(), m))
            table_tex = lowercase(to_latex(m))
            @test occursin("upper scan reached limit", single_plain)
            @test occursin("lower scan reached limit", single_plain)
            @test count("requested confidence crossing: not found", single_plain) == 2
            @test !occursin("best-fit parameter is at", single_plain)
            @test occursin("scan at limit", table_plain)
            @test occursin("boundary displacement", table_plain)
            @test occursin("scan at limit", table_html)
            @test occursin("boundary displacement", table_html)
            @test occursin("\\mathrm{upper\\ scan\\ at\\ limit}", table_tex)
            @test occursin("\\mathrm{lower\\ scan\\ at\\ limit}", table_tex)
            @test occursin("\\mathrm{boundary\\ displacement}", table_tex)

            # Even the compact-value helper must not reinterpret boundary
            # displacements as an ordinary asymmetric MINOS error.
            compact_value = NativeMinuit._value_cell(
                NativeMinuit._param_row_data(m, 1); mode = :text)
            @test occursin("±", compact_value)
            @test !occursin("+", compact_value)
            @test !occursin("−", compact_value)
        end

        @testset "one-sided truncation preserves the crossing side" begin
            m, e = bounded_profile(curvature = 1.0, start = 0.1,
                                   limits = (-0.25, 10.0))
            @test e.upper_valid && e.lower_valid
            @test !e.upper_par_limit && e.lower_par_limit
            @test is_valid(e)
            @test !NativeMinuit.has_closed_interval(e)

            single_plain = sprint(show, MIME"text/plain"(), e)
            table_plain = lowercase(sprint(show, MIME"text/plain"(), m))
            @test occursin("Upper: OK — MINOS error", single_plain)
            @test occursin("Lower scan reached limit", single_plain)
            @test occursin("boundary displacement:", lowercase(single_plain))
            @test occursin("lower scan at limit", table_plain)
            @test occursin("boundary displacement", table_plain)
            @test occursin("+", table_plain)
        end

        @testset "best-fit proximity and scan termination are distinct" begin
            m, e = bounded_profile(
                curvature = 1.0,
                center = 0.05,
                start = 0.1,
                limits = (0.0, 10.0),
            )
            @test 1 in NativeMinuit._at_limit_indices(m)
            @test e.lower_par_limit

            table_plain = sprint(show, MIME"text/plain"(), m)
            table_html = sprint(show, MIME"text/html"(), m)
            @test occursin("[⚠ Parameter(s) close to limit]", table_plain)
            @test occursin("[⚠ MINOS scan reached parameter limit]", table_plain)
            @test occursin("Best-fit value", table_plain)
            @test occursin("MINOS scan reached parameter limit", table_plain)
            @test occursin("Parameter(s) close to limit", table_html)
            @test occursin("MINOS scan reached parameter limit", table_html)
        end
    end

    @testset "closed-interval semantics" begin
        @test is_valid(upper_limit)
        @test is_valid(lower_limit)
        @test is_valid(both_limits)
        @test !NativeMinuit.has_closed_interval(upper_limit)
        @test !NativeMinuit.has_closed_interval(lower_limit)
        @test !NativeMinuit.has_closed_interval(both_limits)
        @test NativeMinuit.has_closed_interval(crossing)
        @test !NativeMinuit.has_closed_interval(upper_invalid)
        @test is_valid(zero_upper)
        @test !NativeMinuit.has_closed_interval(zero_upper)
    end

    @testset "single-result text and LaTeX" begin
        upper_plain = sprint(show, MIME"text/plain"(), upper_limit)
        upper_compact = sprint(show, upper_limit)
        upper_html = sprint(show, MIME"text/html"(), upper_limit)
        @test occursin("upper scan reached limit", lowercase(upper_plain))
        @test occursin("boundary displacement:", lowercase(upper_plain))
        @test occursin("requested confidence crossing: not found", lowercase(upper_plain))
        @test occursin("upper scan at limit", lowercase(upper_compact))
        @test occursin("boundary displacement", lowercase(upper_compact))
        @test occursin("upper scan reached limit", lowercase(upper_html))
        @test occursin("boundary displacement", lowercase(upper_html))

        lower_plain = sprint(show, MIME"text/plain"(), lower_limit)
        both_plain = sprint(show, MIME"text/plain"(), both_limits)
        @test occursin("lower scan reached limit", lowercase(lower_plain))
        @test occursin("boundary displacement", lowercase(lower_plain))
        @test occursin("upper scan reached limit", lowercase(both_plain))
        @test occursin("lower scan reached limit", lowercase(both_plain))

        normal_plain = sprint(show, MIME"text/plain"(), crossing)
        @test !occursin("at limit", lowercase(normal_plain))
        @test occursin("Upper: OK — MINOS error", normal_plain)
        @test occursin("Lower: OK — MINOS error", normal_plain)

        upper_tex = to_latex(upper_limit)
        lower_tex = to_latex(lower_limit)
        both_tex = to_latex(both_limits)
        normal_tex = to_latex(crossing)
        invalid_tex = to_latex(upper_invalid)
        @test occursin("\\mathrm{upper\\ scan\\ at\\ limit}", upper_tex)
        @test occursin("\\mathrm{boundary\\ displacement}", upper_tex)
        @test occursin("\\mathrm{lower\\ scan\\ at\\ limit}", lower_tex)
        @test occursin("\\mathrm{upper\\ scan\\ at\\ limit}", both_tex)
        @test occursin("\\mathrm{lower\\ scan\\ at\\ limit}", both_tex)
        @test !occursin("\\mathrm{scan\\ at\\ limit}", normal_tex)
        @test occursin("^{+", normal_tex) && occursin("_{-", normal_tex)
        @test occursin("\\mathrm{invalid}", invalid_tex)
        @test occursin("_{-", invalid_tex)
    end

    @testset "single-result text box keeps its fixed width" begin
        function box_widths(e)
            text = sprint(show, MIME"text/plain"(), e)
            return textwidth.(split(text, '\n'; keepempty = false))
        end

        cases = (
            "both sides invalid" => both_invalid,
            "upper scan at limit and lower invalid" => upper_limit_lower_invalid,
            "large boundary displacement" => large_upper_limit,
        )
        for (label, e) in cases
            @testset "$label" begin
                widths = box_widths(e)
                @test !isempty(widths)
                @test all(==(73), widths)
            end
        end

        mixed_plain = lowercase(sprint(
            show, MIME"text/plain"(), upper_limit_lower_invalid))
        @test occursin("upper scan reached limit", mixed_plain)
        @test occursin("invalid", mixed_plain)
        @test occursin("boundary displacement", lowercase(sprint(
            show, MIME"text/plain"(), large_upper_limit)))
    end

    @testset "fit-table text, HTML, and LaTeX" begin
        target = min_value
        m = Minuit(
            p -> (p[1] - target)^2 + (p[2] - 1.0)^2,
            [3.7, 0.0];
            names = ["x", "y"],
            limits = [(nothing, 3.78), nothing],
        )
        migrad!(m)
        m.minos_errors[1] = MinosError(
            1, m.values[1], 3.78 - m.values[1], -0.3,
            true, true, false, false, false, false,
            true, false, 17,
        )

        table_plain = sprint(show, MIME"text/plain"(), m)
        table_html = sprint(show, MIME"text/html"(), m)
        table_tex = to_latex(m)
        @test occursin("upper scan at limit", lowercase(table_plain))
        @test occursin("boundary displacement", lowercase(table_plain))
        @test occursin("upper scan at limit", lowercase(table_html))
        @test occursin("boundary displacement", lowercase(table_html))
        @test occursin("\\mathrm{upper\\ scan\\ at\\ limit}", table_tex)
        @test occursin("\\mathrm{boundary\\ displacement}", table_tex)
        @test occursin("MINOS scan reached parameter limit", table_plain)
        @test occursin("MINOS scan reached parameter limit", table_html)
    end

    @testset "plot recipes omit non-crossing whiskers" begin
        function recipe_yerror(e)
            rd = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(), e))
            return rd.plotattributes[:yerror], rd.plotattributes[:label]
        end

        yerr_upper, label_upper = recipe_yerror(upper_limit)
        @test yerr_upper == ([0.3], [0.0])
        @test occursin("interval truncated by parameter limit", lowercase(label_upper))
        @test occursin("upper boundary displacement omitted", lowercase(label_upper))

        yerr_lower, label_lower = recipe_yerror(lower_limit)
        @test yerr_lower == ([0.0], [0.2])
        @test occursin("interval truncated by parameter limit", lowercase(label_lower))
        @test occursin("lower boundary displacement omitted", lowercase(label_lower))

        yerr_both, label_both = recipe_yerror(both_limits)
        @test yerr_both == ([0.0], [0.0])
        @test occursin("interval truncated by parameter limits", lowercase(label_both))
        @test occursin("boundary displacements omitted", lowercase(label_both))

        yerr_crossing, _ = recipe_yerror(crossing)
        @test yerr_crossing == ([0.3], [0.2])

        yerr_invalid, label_invalid = recipe_yerror(upper_invalid)
        @test yerr_invalid == ([0.3], [0.0])
        @test occursin("invalid side omitted", lowercase(label_invalid))

        vector_recipe = only(RecipesBase.apply_recipe(
            Dict{Symbol,Any}(), [upper_limit, lower_limit, both_limits, crossing, upper_invalid]))
        @test vector_recipe.plotattributes[:yerror] ==
              ([0.3, 0.0, 0.0, 0.3, 0.3], [0.0, 0.2, 0.0, 0.2, 0.0])
        @test occursin("intervals truncated by parameter limits",
                       lowercase(vector_recipe.plotattributes[:label]))
    end

    @testset "serialization preserves boundary flags and compatibility" begin
        for e in (upper_limit, lower_limit, both_limits, crossing, upper_invalid)
            d = NativeMinuit.to_dict(e)
            @test d["upper"] == e.upper
            @test d["lower"] == e.lower
            @test d["upper_valid"] == e.upper_valid
            @test d["lower_valid"] == e.lower_valid
            @test d["upper_par_limit"] == e.upper_par_limit
            @test d["lower_par_limit"] == e.lower_par_limit

            restored = NativeMinuit.minos_error_from_dict(JSON.parse(JSON.json(d)))
            @test restored.upper == e.upper
            @test restored.lower == e.lower
            @test restored.upper_valid == e.upper_valid
            @test restored.lower_valid == e.lower_valid
            @test restored.upper_par_limit == e.upper_par_limit
            @test restored.lower_par_limit == e.lower_par_limit
            @test NativeMinuit.has_closed_interval(restored) ==
                  NativeMinuit.has_closed_interval(e)
        end

        legacy = NativeMinuit.to_dict(crossing)
        delete!(legacy, "upper_par_limit")
        delete!(legacy, "lower_par_limit")
        restored_legacy = NativeMinuit.minos_error_from_dict(legacy)
        @test !restored_legacy.upper_par_limit
        @test !restored_legacy.lower_par_limit
    end
end
