module SimsFlanagan

using LinearAlgebra
using StaticArrays
using Random

# SciMLBase interface - we extend these functions
using SciMLBase
import SciMLBase: solve, remake

# Astrodynamics
using AstroCoords
using Lambert: Lambert  # import to avoid solve conflict

# Optimization - import to avoid solve conflict  
using Optimization: Optimization
using OptimizationMOI
using OptimizationMadNLP
using MadNLP: MadNLP
using ForwardDiff

include("utils.jl")
include("types.jl")
include("propagation.jl")
include("problem.jl")
include("solve.jl")

end # module SimsFlanagan
