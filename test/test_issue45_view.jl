# SPDX-License-Identifier: LGPL-2.1-or-later

# Issue #45 (design §6.2) — the `p.pars` derived read-only view, the
# `values=` donation semantics of the record constructor, the exception
# surface of every mutation entry point, and generic-container acceptance
# beyond ComponentArrays.
#
# These tests are the AUTHORITY behind the CHANGELOG/docstring wording on
# exception types (design §4; round-2 glm #R2-1): `setindex!` and `sort!`
# with a `by=`/`lt=` ordering throw the guiding `ArgumentError`, while bare
# `sort!` (no `isless` for `MinuitParameter`) and the growing/shrinking
# mutators (`push!`/`resize!`/`empty!`/`append!`) throw `MethodError`.

# ── top-level helpers (included files evaluate at module top level) ──────────

# Function barrier for the allocation gate: sums `p.pars[i].value` the way a
# hot loop would. Defined at top level so `@allocated` measures the call
# itself, not testset-scope capture artifacts (recorded v0.3.0 trap:
# `@allocated` at testset scope needs function barriers).
function _issue45_sum_pars_values(p::Parameters)
    s = 0.0
    for i in 1:n_pars(p)
        s += p.pars[i].value
    end
    return s
end
_issue45_alloc_sum(p::Parameters) = @allocated _issue45_sum_pars_values(p)

# Minimal custom container (design §6.2 "generic-container acceptance"): an
# `AbstractVector{Float64}` that is NOT a `Vector` and whose `similar` (the
# `Base` fallback) returns a plain `Vector{Float64}`. `Parameters` must bind
# its container type parameter from the ACTUAL `similar` allocation, not
# from `typeof(values)`.
struct Issue45LocalVec <: AbstractVector{Float64}
    data::Vector{Float64}
end
Base.size(v::Issue45LocalVec) = size(v.data)
Base.IndexStyle(::Type{Issue45LocalVec}) = IndexLinear()
Base.getindex(v::Issue45LocalVec, i::Int) = v.data[i]
Base.setindex!(v::Issue45LocalVec, x, i::Int) = (v.data[i] = x; v)

@testset "issue #45: p.pars view + donation + exception surface" begin

    # Bounded + fixed configuration exercised throughout: a two-sided bound,
    # an unbounded free parameter, and a lower-bounded FIXED parameter.
    mkpars() = [
        MinuitParameter("alpha", 1.5, 0.1; lower = 0.0, upper = 4.0),
        MinuitParameter("beta", -0.5, 0.2),
        MinuitParameter("gamma", 2.5, 0.3; lower = 1.0, fixed = true),
    ]

    @testset "view semantics: eltype/length/collect/iteration/roundtrip" begin
        pars = mkpars()
        P = Parameters(pars)
        v = P.pars

        @test v isa AbstractVector{MinuitParameter}
        @test eltype(v) == MinuitParameter
        @test length(v) == 3
        @test size(v) == (3,)
        @test collect(v) isa Vector{MinuitParameter}
        @test length(collect(v)) == 3

        # Iteration materializes the records in order.
        @test [p.name for p in v] == ["alpha", "beta", "gamma"]

        # Full roundtrip of every record field against the canonical
        # stores: numbers from `values`/`errors`, structure from metadata
        # (compared here against the constructor-time records).
        for i in 1:3
            r = v[i]
            @test r isa MinuitParameter
            @test r.name == pars[i].name
            @test r.value == P.values[i]
            @test r.error == P.errors[i]
            @test r.value == pars[i].value
            @test r.error == pars[i].error
            @test isequal(r.lower, pars[i].lower)   # NaN-safe
            @test isequal(r.upper, pars[i].upper)
            @test r.fixed == pars[i].fixed
        end
        # The derived record is egal to the source record (immutable
        # struct: `===` compares field-wise, NaN bitwise).
        @test v[3] === pars[3]

        # The view is derived, not stored: it reflects the canonical
        # vectors at READ time.
        getfield(P, :values)[2] = -0.25
        @test P.pars[2].value == -0.25
    end

    @testset "inference, zero allocation, ===-egal reads" begin
        P = Parameters(mkpars())

        @test (@inferred P.pars[1]) isa MinuitParameter

        # Warm up, then pin: the `p.pars[i].value` idiom behind a function
        # barrier must not allocate (the `@inline getproperty` +
        # immutable-wrapper design, §3.3/§5).
        _issue45_sum_pars_values(P)
        _issue45_alloc_sum(P)
        @test _issue45_alloc_sum(P) == 0
        @test _issue45_sum_pars_values(P) == 1.5 + (-0.5) + 2.5

        # Consecutive reads of the same index are egal — the record is a
        # pure function of the canonical stores.
        @test P.pars[1] === P.pars[1]
        @test P.pars[3] === P.pars[3]
    end

    @testset "exception surface pins (design §4 authority)" begin
        P = Parameters(mkpars())
        v = P.pars
        rec = MinuitParameter("zz", 0.0, 1.0)

        # setindex! → guiding ArgumentError naming the supported mutators.
        err = try
            v[1] = rec
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("set_value!", err.msg)

        # sort! WITH an ordering (comparisons succeed on Float64, first
        # in-place write routes through setindex!) → same ArgumentError.
        # Values are deliberately out of ascending order so the sort cannot
        # early-exit before attempting a write.
        @test_throws ArgumentError sort!(v; by = p -> p.value)

        # BARE sort!: no `isless(::MinuitParameter, ::MinuitParameter)` —
        # the comparison itself is the missing method.
        @test_throws MethodError sort!(v)

        # Growing/shrinking mutators: no resizing methods exist for the view.
        @test_throws MethodError push!(v, rec)
        @test_throws MethodError resize!(v, 5)
        @test_throws MethodError empty!(v)
        @test_throws MethodError append!(v, [rec])

        # None of the failed attempts touched the canonical stores.
        @test [p.name for p in P.pars] == ["alpha", "beta", "gamma"]
        @test collect(P.values) == [1.5, -0.5, 2.5]
        @test collect(P.errors) == [0.1, 0.2, 0.3]
    end

    @testset "curated propertynames" begin
        P = Parameters(mkpars())
        @test propertynames(P) == (:pars, :values, :errors, :prec)
        @test :metadata in propertynames(P, true)
        @test :pars in propertynames(P, true)
        # Internal fields stay reachable but are not advertised publicly.
        @test :metadata ∉ propertynames(P)
    end

    @testset "donation semantics: record numbers always win" begin
        # `values=` donates ONLY the container/axes — the numbers for BOTH
        # canonical vectors come from the records (round-1 glm #8a).
        P = Parameters([MinuitParameter("x", 1.0, 0.5)], MachinePrecision();
                       values = [99.0])
        @test collect(P.values) == [1.0]
        @test collect(P.errors) == [0.5]

        # `Parameters(pars, source)` takes the errors from `pars[i].error`,
        # NOT from `source.errors` — pins the `_build_resume_params`
        # floored-resume-step semantics (design §3.2).
        source = Parameters([MinuitParameter("x", 7.0, 9.0),
                             MinuitParameter("y", 8.0, 9.0)])
        recs = [MinuitParameter("x", 1.0, 0.5),
                MinuitParameter("y", 2.0, 0.25)]
        P2 = Parameters(recs, source)
        @test collect(P2.values) == [1.0, 2.0]     # records, not source values
        @test collect(P2.errors) == [0.5, 0.25]    # records, not source errors
        @test P2.prec === source.prec
        # The donor itself is untouched.
        @test collect(source.values) == [7.0, 8.0]
        @test collect(source.errors) == [9.0, 9.0]
    end

    @testset "generic container beyond ComponentArrays" begin
        recs = [MinuitParameter("u", 1.0, 0.1),
                MinuitParameter("v", 2.0, 0.2)]
        donor = Issue45LocalVec([50.0, 60.0])
        Plv = Parameters(recs, MachinePrecision(); values = donor)

        # Type parameter bound from the ACTUAL allocation: `similar` on the
        # custom container returns a plain Vector, so P === Vector{Float64}.
        @test Plv isa Parameters{Vector{Float64}}
        @test getfield(Plv, :values) isa Vector{Float64}
        @test getfield(Plv, :errors) isa Vector{Float64}
        # Numbers from the records; the donor keeps its own numbers.
        @test Plv.values == [1.0, 2.0]
        @test Plv.errors == [0.1, 0.2]
        @test collect(donor) == [50.0, 60.0]
        # The view works over the generic-container build too.
        @test Plv.pars[2].value == 2.0
        @test Plv.pars[2].error == 0.2

        # End-to-end sanity: a tiny migrad! fit on plain vectors still runs
        # and converges through the new storage layout.
        m = Minuit(x -> (x[1] - 1.0)^2 + (x[2] - 2.0)^2, [0.0, 0.0];
                   errors = [0.1, 0.1])
        migrad!(m)
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-3
        @test m.values[2] ≈ 2.0 atol = 1e-3
        @test m.params.pars[1].value == m.values[1]
    end
end
