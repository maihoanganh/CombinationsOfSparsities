module CombinationsOfSparsities

using MultivariatePolynomials
using DynamicPolynomials
using JuMP
using MosekTools
using SparseArrays
using LinearAlgebra
using LightGraphs
using RowEchelon

export extract_optimizers, Pu_block_hierarchy, Pu_sparse_hierarchy, Pu_mix_hierarchy, PuVa_block_hierarchy, PuVa_sparse_hierarchy, PuVa_mix_hierarchy, block_ASC, sparse_ASC, mix_ASC

include("based_functions.jl")
include("term_sparsity.jl")
include("correlative_sparsity.jl")
include("mix_term_correlative_sparsity.jl")

end # module
