module EllipsotopeModule

using Reexport, Requires
using Random: AbstractRNG, GLOBAL_RNG
using ReachabilityBase.Distribution: reseed!
using ReachabilityBase.Require: require
using LinearAlgebra: norm, dot

using ..LazySets: AbstractCentrallySymmetric, LazySets, AbstractZonotope, Zonotope, Ellipsoid,
                ngens, genmat, norm, constraints_list
using ..LazySets.JuMP: Model, @variable, @objective, @constraint, optimize!, set_silent,
              termination_status, objective_value, value, SecondOrderCone, OPTIMAL, ALMOST_OPTIMAL

@reexport import ..LazySets: dim, rand, genmat, ngens, norm, constraints_list, in
@reexport import ..API: center, ∈, isempty, ρ

export Ellipsotope

include("Ellipsotope.jl")
include("genmat.jl")
include("ngens.jl")
include("center.jl")
include("in.jl")
include("dim.jl")
include("norm.jl")
include("constraints_list.jl")
include("isempty.jl")
include("support_function.jl")

include("init.jl")

end