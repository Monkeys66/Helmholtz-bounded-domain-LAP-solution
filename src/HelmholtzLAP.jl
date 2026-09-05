module HelmholtzLAP

using Gridap
using GridapGmsh
using LinearAlgebra
using Arpack
using IterativeSolvers
using Printf

include("fem.jl")
export assemble_fem
export assemble_source_loads

include("resonance.jl")
export compute_eigenpairs
export show_eigenvalues
export find_resonant_indices
export check_source_orthogonality
export project_load_to_range

include("LAP.jl")
export prepare_source_extension
export compute_particular_solution
export compute_lap_solution

end
