module CombinationsOfSparsities

using MultivariatePolynomials
using DynamicPolynomials
using JuMP
using MosekTools
using LinearAlgebra
using LightGraphs
using RowEchelon

export get_basis, bfind, get_blocks, info, extract_optimizers, Pu_block_hierarchy, Pu_sparse_hierarchy, Pu_mix_hierarchy

include("based_functions.jl")
include("term_sparsity.jl")
include("correlative_sparsity.jl")
include("mix_term_correlative_sparsity.jl")

end # module
