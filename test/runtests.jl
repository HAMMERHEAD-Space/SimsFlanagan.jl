using Test
using SimsFlanagan
using LinearAlgebra
using Random
using StaticArrays

@testset "SimsFlanagan.jl Tests" begin
    include("test_types.jl")
    include("test_propagation.jl")
    include("test_problem.jl")
    include("test_solve.jl")
    include("test_propulsion_types.jl")
    include("test_sundman.jl")
    include("test_pykep_validation.jl")
end
