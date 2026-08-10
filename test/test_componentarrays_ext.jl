# SPDX-License-Identifier: LGPL-2.1-or-later

using ComponentArrays

@testset "ComponentArrays extension" begin
    start = ComponentArray(a = 0.0, b = 0.0)
    expected_axes = ComponentArrays.getaxes(start)

    @testset "objective receives the original axes" begin
        saw_axes = Ref(false)
        borrowed_parent = Ref{Any}()
        objective(p) = begin
            saw_axes[] = ComponentArrays.getaxes(p) == expected_axes
            borrowed_parent[] = parent(p)
            return (p.a - 1.0)^2 + (p.b - 2.0)^2
        end

        m = Minuit(objective, start; errors = [0.1, 0.1])

        # The labeled callback view borrows the minimizer workspace rather
        # than copying it.
        probe = [3.0, 4.0]
        @test @inferred(m.fcn.f(probe)) == 8.0
        @test borrowed_parent[] === probe

        migrad!(m)
        @test saw_axes[]
        @test m.valid
        @test m.values[1] ≈ 1.0 atol = 1e-4
        @test m.values[2] ≈ 2.0 atol = 1e-4
        @test m.fmin.ext_values isa Vector{Float64}
        @test m.values isa ComponentVector
        @test m.values.a ≈ 1.0 atol = 1e-4
        @test m.values.b ≈ 2.0 atol = 1e-4
        @test ComponentArrays.getaxes(m.values) == expected_axes
        @test ComponentArrays.getaxes(m.errors) == expected_axes

        # Labeled result views retain NativeMinuit's write-through semantics,
        # including fit-cache invalidation, and copy construction keeps the
        # original external representation.
        m.values.a = 0.25
        @test m.fmin === nothing
        @test m.values.a == 0.25
        m_copy = Minuit(objective, m)
        @test m_copy.values isa ComponentVector
        @test m_copy.values.a == 0.25
        @test ComponentArrays.getaxes(m_copy.values) == expected_axes

        # Bound transforms still reconstruct the full external vector before
        # adapting it for the user callback.
        mb = Minuit(objective, start;
                    errors = [0.1, 0.1], limits = [(-5.0, 5.0), nothing])
        migrad!(mb)
        @test mb.valid
        @test mb.values[1] ≈ 1.0 atol = 1e-2
        @test mb.values[2] ≈ 2.0 atol = 1e-4
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
