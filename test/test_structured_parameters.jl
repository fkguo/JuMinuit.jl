# SPDX-License-Identifier: LGPL-2.1-or-later

using ComponentArrays

@testset "structured parameter containers" begin
    start = ComponentArray(a = 0.0, b = 0.0)
    expected_axes = ComponentArrays.getaxes(start)

    @testset "objective receives the original axes" begin
        saw_axes = Ref(false)
        objective(p) = begin
            saw_axes[] = ComponentArrays.getaxes(p) == expected_axes
            return (p.a - 1.0)^2 + (p.b - 2.0)^2
        end

        m = Minuit(objective, start; errors = [0.1, 0.1])

        # No dependency-specific callback wrapper: the user's function is
        # retained directly, while full external workspaces follow x0's
        # ordinary `similar` behavior.
        @test m.fcn.f === objective
        @test m.params.prototype isa ComponentVector
        @test NativeMinuit.int_to_ext_vector(
            m.params, NativeMinuit.initial_int_values(m.params)) isa ComponentVector

        migrad!(m)
        @test saw_axes[]
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-4
        @test m.values[2] ≈ 2.0 atol = 1e-4
        @test m.fmin.ext_values isa ComponentVector
        @test m.fmin.ext_errors isa ComponentVector
        @test m.fmin.ext_values.a ≈ 1.0 atol = 1e-4
        @test m.fmin.ext_values.b ≈ 2.0 atol = 1e-4
        @test ComponentArrays.getaxes(m.fmin.ext_values) == expected_axes
        @test ComponentArrays.getaxes(m.fmin.ext_errors) == expected_axes

        # Copy construction keeps the external container through the same
        # generic prototype path.
        m_copy = Minuit(objective, m)
        migrad!(m_copy)
        @test m_copy.fmin.ext_values isa ComponentVector
        @test ComponentArrays.getaxes(m_copy.fmin.ext_values) == expected_axes

        # Parameter mutation rebuilds metadata without discarding the external
        # container prototype.
        set_value!(m_copy, 1, 0.5)
        @test m_copy.params.prototype isa ComponentVector
        migrad!(m_copy)
        @test m_copy.valid

        # Bound transforms still reconstruct the full structured external
        # vector before calling the user objective.
        mb = Minuit(objective, start;
                    errors = [0.1, 0.1], limits = [(-5.0, 5.0), nothing])
        migrad!(mb)
        @test mb.valid
        @test mb.values[1] ≈ 1.0 atol = 1e-2
        @test mb.values[2] ≈ 2.0 atol = 1e-4
    end

    @testset "nested axes are generic" begin
        nested_start = ComponentArray(
            line = (slope = 0.0, intercept = 0.0),
            offset = 0.0,
        )
        objective(p) = (p.line.slope - 1.0)^2 +
                       (p.line.intercept - 2.0)^2 +
                       (p.offset - 3.0)^2

        m = Minuit(objective, nested_start; errors = fill(0.1, 3))
        migrad!(m)

        @test m.valid
        @test m.fmin.ext_values.line.slope ≈ 1.0 atol = 1e-4
        @test m.fmin.ext_values.line.intercept ≈ 2.0 atol = 1e-4
        @test m.fmin.ext_values.offset ≈ 3.0 atol = 1e-4
    end

    @testset "user gradient receives the original axes" begin
        saw_objective_axes = Ref(false)
        saw_gradient_axes = Ref(false)
        objective(p) = begin
            saw_objective_axes[] = ComponentArrays.getaxes(p) == expected_axes
            return (p.a - 1.0)^2 + (p.b - 2.0)^2
        end
        gradient(p) = begin
            saw_gradient_axes[] = ComponentArrays.getaxes(p) == expected_axes
            return ComponentArray(a = 2 * (p.a - 1.0),
                                  b = 2 * (p.b - 2.0))
        end

        m = Minuit(objective, start;
                   errors = [0.1, 0.1], grad = gradient)
        migrad!(m)

        @test saw_objective_axes[]
        @test saw_gradient_axes[]
        @test m.ngrad > 0
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-4
        @test m.values[2] ≈ 2.0 atol = 1e-4
    end
end
