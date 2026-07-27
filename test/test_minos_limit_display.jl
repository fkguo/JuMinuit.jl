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

    @testset "closed-interval semantics" begin
        @test is_valid(upper_limit)
        @test is_valid(lower_limit)
        @test is_valid(both_limits)
        @test !NativeMinuit.has_closed_interval(upper_limit)
        @test !NativeMinuit.has_closed_interval(lower_limit)
        @test !NativeMinuit.has_closed_interval(both_limits)
        @test NativeMinuit.has_closed_interval(crossing)
        @test !NativeMinuit.has_closed_interval(upper_invalid)
    end

    @testset "single-result text and LaTeX" begin
        upper_plain = sprint(show, MIME"text/plain"(), upper_limit)
        upper_compact = sprint(show, upper_limit)
        @test occursin("upper at limit", lowercase(upper_plain))
        @test occursin("distance", lowercase(upper_plain))
        @test occursin("upper at limit", lowercase(upper_compact))
        @test occursin("distance", lowercase(upper_compact))

        lower_plain = sprint(show, MIME"text/plain"(), lower_limit)
        both_plain = sprint(show, MIME"text/plain"(), both_limits)
        @test occursin("lower at limit", lowercase(lower_plain))
        @test occursin("distance", lowercase(lower_plain))
        @test occursin("upper at limit", lowercase(both_plain))
        @test occursin("lower at limit", lowercase(both_plain))

        normal_plain = sprint(show, MIME"text/plain"(), crossing)
        @test !occursin("at limit", lowercase(normal_plain))
        @test occursin("Upper: OK", normal_plain)
        @test occursin("Lower: OK", normal_plain)

        upper_tex = to_latex(upper_limit)
        lower_tex = to_latex(lower_limit)
        both_tex = to_latex(both_limits)
        normal_tex = to_latex(crossing)
        invalid_tex = to_latex(upper_invalid)
        @test occursin("\\mathrm{at\\ upper\\ limit}", upper_tex)
        @test occursin("\\mathrm{distance}", upper_tex)
        @test occursin("\\mathrm{at\\ lower\\ limit}", lower_tex)
        @test occursin("\\mathrm{at\\ upper\\ limit}", both_tex)
        @test occursin("\\mathrm{at\\ lower\\ limit}", both_tex)
        @test !occursin("\\mathrm{at\\ ", normal_tex)
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
            "upper at limit and lower invalid" => upper_limit_lower_invalid,
            "large at-limit distance" => large_upper_limit,
        )
        for (label, e) in cases
            @testset "$label" begin
                widths = box_widths(e)
                @test length(widths) == 9
                @test all(==(73), widths)
            end
        end

        mixed_plain = lowercase(sprint(
            show, MIME"text/plain"(), upper_limit_lower_invalid))
        @test occursin("at upper limit", mixed_plain)
        @test occursin("invalid", mixed_plain)
        @test occursin("distance", lowercase(sprint(
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
        @test occursin("upper at limit", lowercase(table_plain))
        @test occursin("distance", lowercase(table_plain))
        @test occursin("upper at limit", lowercase(table_html))
        @test occursin("distance", lowercase(table_html))
        @test occursin("\\mathrm{at\\ upper\\ limit}", table_tex)
        @test occursin("\\mathrm{distance}", table_tex)
    end

    @testset "plot recipes omit non-crossing whiskers" begin
        function recipe_yerror(e)
            rd = only(RecipesBase.apply_recipe(Dict{Symbol,Any}(), e))
            return rd.plotattributes[:yerror], rd.plotattributes[:label]
        end

        yerr_upper, label_upper = recipe_yerror(upper_limit)
        @test yerr_upper == ([0.3], [0.0])
        @test occursin("upper at limit", lowercase(label_upper))

        yerr_lower, label_lower = recipe_yerror(lower_limit)
        @test yerr_lower == ([0.0], [0.2])
        @test occursin("lower at limit", lowercase(label_lower))

        yerr_both, label_both = recipe_yerror(both_limits)
        @test yerr_both == ([0.0], [0.0])
        @test occursin("limits omitted", lowercase(label_both))

        yerr_crossing, _ = recipe_yerror(crossing)
        @test yerr_crossing == ([0.3], [0.2])

        yerr_invalid, label_invalid = recipe_yerror(upper_invalid)
        @test yerr_invalid == ([0.3], [0.0])
        @test occursin("invalid side omitted", lowercase(label_invalid))

        vector_recipe = only(RecipesBase.apply_recipe(
            Dict{Symbol,Any}(), [upper_limit, lower_limit, both_limits, crossing, upper_invalid]))
        @test vector_recipe.plotattributes[:yerror] ==
              ([0.3, 0.0, 0.0, 0.3, 0.3], [0.0, 0.2, 0.0, 0.2, 0.0])
        @test occursin("non-crossing sides omitted",
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
