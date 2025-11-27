module SimsFlanagan

using LinearAlgebra
using StaticArrays

using AstroCoords
using Lambert

using Optimization
using OptimizationOptimJL

include("types.jl")
include("propagation.jl")
include("problem.jl")
include("solve.jl")

end # module SimsFlanagan
