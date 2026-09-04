module HelmholtzLAP

using Gridap
using GridapGmsh
using LinearAlgebra
using Arpack
using IterativeSolvers
using Printf

include("fem.jl")

export assemble_fem
export compute_eigenpairs
export show_eigens
export compute_parSolution

end

