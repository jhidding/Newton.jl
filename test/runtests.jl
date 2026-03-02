module Spec

using Newton, Test, Unitful
using Newton: set_still!, random_particles, momentum, MomentumVector, run_simulation,
    total_energy

@testset "Newton.jl" begin
    @testset "testing sanity" begin
        @test 1 + 1 == 2
    end

    @testset "run $i" for i in 1:10
        p = set_still!(random_particles(3, seed=i))
        @test momentum(p) ≈ zero(MomentumVector) atol=1e-6u"N*s"

        orbit = run_simulation(p, 0.1u"s", 1000)
        @test momentum(orbit[end]) ≈ zero(MomentumVector) atol=1e-6u"N*s"
        @test total_energy(orbit[1]) ≈ total_energy(orbit[end]) rtol=1e-6
    end
end

end  # module Spec
