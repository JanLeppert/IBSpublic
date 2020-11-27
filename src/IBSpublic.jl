module IBSpublic
# using the module developed in github "TemperatureGradientGCibs.jl" 
# from 29.10.2020 (Branch 'dimensionless', commit 'include an old version of Solving.jl', fc4fc4c)

__precompile__()
using StaticArrays
using DifferentialEquations
include("Types.jl")
using Dierckx
include("Model.jl")
using ForwardDiff
using Interpolations
include("Solving.jl")

end
