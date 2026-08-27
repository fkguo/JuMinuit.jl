# SPDX-License-Identifier: LGPL-2.1-or-later

# Issue #45 (design §6.2) — fit-overlay guarantee, overlay isolation over
# ALL SIX copied containers (metadata, ext_of_int, int_of_ext, name_to_ext,
# values, errors), post-`reset` exception atomicity of bulk setters, and the
# issue-#46 fix!/release! chain through the new metadata/values storage.

@testset "issue #45: overlay isolation + atomicity" begin

    # Bounded + fixed three-parameter fit used throughout: "a" two-sided
    # bounded, "b" unbounded, "c" FIXED (at its own optimum, so the fit is
    # clean). Minimum at (1, 2, 3).
    mk3() = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2 + (x[3] - 3.0)^2,
                   [0.5, 0.0, 3.0];
                   names = ["a", "b", "c"], errors = [0.1, 0.1, 0.1],
                   limits = [(0.0, 5.0), nothing, nothing],
                   fixed = [false, false, true])

    @testset "fit-overlay guarantee: pars view tracks values/errors" begin
        m = mk3()

        # Before any fit: the config values/steps.
        for i in 1:3
            @test m.params.pars[i].value == m.values[i]
            @test m.params.pars[i].error == m.errors[i]
        end
        @test m.values[3] == 3.0

        # After migrad!: the overlay reflects the FIT.
        migrad!(m)
        @test m.valid
        for i in 1:3
            @test m.params.pars[i].value == m.values[i]
            @test m.params.pars[i].error == m.errors[i]
        end
        # Bounded "a" stops on EDM under the sin-transform — 1e-2 like the
        # bounded case in test_minuit_mutators.jl.
        @test m.values[1] ≈ 1.0 atol = 1e-2
        @test m.values[2] ≈ 2.0 atol = 1e-3
        @test m.values[3] == 3.0                      # fixed stays put
        @test m.params.pars[3].fixed

        # After minos! (runs on the free parameters): overlay still in sync.
        minos!(m)
        @test length(m.minos_errors) == 2
        for i in 1:3
            @test m.params.pars[i].value == m.values[i]
            @test m.params.pars[i].error == m.errors[i]
        end
    end

    @testset "overlay isolation over all six copied containers" begin
        m = mk3()
        migrad!(m)
        @test m.valid

        # Baselines captured BEFORE the sabotage.
        fit_vals = collect(m.values)
        fit_errs = collect(m.errors)
        cfg = NativeMinuit._init_params(m)
        cfg_vals = copy(getfield(cfg, :values))
        cfg_errs = copy(getfield(cfg, :errors))
        fmin_vals = copy(m.fmin.ext_values)

        # Repeated m.params calls return independent objects post-fit.
        @test m.params !== m.params

        # Mutate EVERY mutable container reachable from one overlay
        # (§3.4 invariant: omitting any single `copy(...)` in
        # `_fit_overlaid_params` must fail at least one assertion below).
        ov = m.params
        ov.values[1] = -777.0                          # 1: values
        ov.errors[1] = -777.0                          # 2: errors
        reverse!(getfield(ov, :metadata))              # 3: metadata
        ov.name_to_ext["zzz_sabotage"] = 99            # 4: name_to_ext
        ov.ext_of_int[1] = 42                          # 5: ext_of_int
        ov.int_of_ext[1] = 42                          # 6: int_of_ext

        # A FRESH m.params is unaffected in every container.
        fresh = m.params
        @test fresh !== ov
        @test collect(fresh.values) == fit_vals
        @test collect(fresh.errors) == fit_errs
        @test [md.name for md in getfield(fresh, :metadata)] == ["a", "b", "c"]
        @test !haskey(fresh.name_to_ext, "zzz_sabotage")
        @test getfield(fresh, :ext_of_int) == [1, 2]   # "c" is fixed
        @test getfield(fresh, :int_of_ext) == [1, 2, 0]
        for i in 1:3
            @test fresh.pars[i].value == m.values[i]
        end

        # The stored config internals are unaffected.
        cfg2 = NativeMinuit._init_params(m)
        @test collect(getfield(cfg2, :values)) == cfg_vals
        @test collect(getfield(cfg2, :errors)) == cfg_errs
        @test [md.name for md in getfield(cfg2, :metadata)] == ["a", "b", "c"]
        @test !haskey(getfield(cfg2, :name_to_ext), "zzz_sabotage")
        @test getfield(cfg2, :ext_of_int) == [1, 2]
        @test getfield(cfg2, :int_of_ext) == [1, 2, 0]

        # The cached fit is unaffected.
        @test collect(m.fmin.ext_values) == fmin_vals

        # A subsequent migrad! (resume from the cached fit) is unaffected.
        migrad!(m)
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-2    # bounded: EDM-stop tolerance
        @test m.values[2] ≈ 2.0 atol = 1e-3
        @test m.values[3] == 3.0
    end

    @testset "post-reset exception atomicity of a failing bulk set" begin
        m = mk3()
        # Commit a real config edit first, so "the pre-attempt committed
        # config" is distinguishable from the constructor state.
        set_value!(m, 2, 0.25)
        migrad!(m)
        @test m.valid

        # Bulk set with one invalid entry: throws, commits NOTHING (the
        # per-element validation shares `_build_value_par` with the
        # per-parameter mutator; commit is a single `_replace_all_params!`).
        @test_throws ArgumentError (m.values = [1.0, NaN, 3.0])
        @test m.fmin !== nothing                     # fit not dropped either

        # As seen from a fresh cold start: the committed config — the
        # constructor state plus the ONE legitimate edit — not a partial
        # application of the failed bulk assignment.
        reset(m)
        @test m.fmin === nothing
        @test collect(m.values) == [0.5, 0.25, 3.0]
        @test collect(m.errors) == [0.1, 0.1, 0.1]
        @test isequal(collect(m.limits),
                      [(0.0, 5.0), (NaN, NaN), (NaN, NaN)])
        @test collect(m.fixed) == [false, false, true]
        # Maps unchanged too (would differ had fixed flags been touched).
        @test n_free(NativeMinuit._init_params(m)) == 2
        @test getfield(NativeMinuit._init_params(m), :int_of_ext) == [1, 2, 0]
    end

    @testset "issue #46 chain through the new storage" begin
        m = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2, [0.0, 0.0];
                   names = ["a", "b"], errors = [0.1, 0.1])
        migrad!(m)
        @test m.valid
        v_fit = collect(m.values)
        @test v_fit[1] ≈ 1.0 atol = 1e-3
        @test v_fit[2] ≈ 2.0 atol = 1e-3

        # fix! between fits freezes par 1 AT ITS FITTED value; par 2 keeps
        # its fitted value; the steps stay at the user-set ones (values-only
        # commit — retry length scale / resume floor preserved).
        fix!(m, 1)
        @test m.fmin === nothing
        @test m.values[1] == v_fit[1]
        @test m.values[2] == v_fit[2]
        @test collect(m.errors) == [0.1, 0.1]
        @test m.fixed[1] == true

        release!(m, 1)
        @test m.fixed[1] == false
        @test collect(m.values) == v_fit

        # The chain ends in a working fit from the current point.
        migrad!(m)
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-3
        @test m.values[2] ≈ 2.0 atol = 1e-3
    end
end
