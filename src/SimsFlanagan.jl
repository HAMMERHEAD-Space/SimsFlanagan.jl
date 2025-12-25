module SimsFlanagan

using LinearAlgebra
using StaticArrays
using Random

# SciMLBase interface - we extend these functions
using SciMLBase
import SciMLBase: solve, remake

# Astrodynamics
using AstroCoords
import Lambert  # import to avoid solve conflict

# Optimization - import to avoid solve conflict  
import Optimization
using OptimizationMOI
using OptimizationMadNLP
import MadNLP
import MadNLPMumps
using ForwardDiff

include("utils.jl")
include("types.jl")
include("propagation.jl")
include("problem.jl")
include("solve.jl")

end # module SimsFlanagan
