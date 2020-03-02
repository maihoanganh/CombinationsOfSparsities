include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    eps=1e-5
    k=2 # relaxed order
    r=10 # sparse order

    @polyvar x1 x2 x3 # define polnomial variables
    x=[x1;x2;x3]

    f = x1+x2+x3 # objective function
    g = [x1;x2;x3] # inequality constraints
    h = [x1*x2*x3-1] # equality constraints

    opt_val,sol = CombinationsOfSparsities.PuVa_block_hierarchy(x,f,g,h,eps,k,r,knownlb=0)

    g=[g;opt_val-f]

    eps=1e-2
    k=1 # relaxed order
    r=2 # sparse order

    sol=CombinationsOfSparsities.term_sparsity_adding_spherical_constraints(x,g,h,eps,k,r) # compute minimizers
end


"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]]]
Termination status = OPTIMAL
Optimal value = 3.0000398046112338
Rank of moment matrix = 11
Determine L0:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001]]]
Termination status = OPTIMAL
Optimal value = 3.039999891198467
Rank of moment matrix = 1
------------------------------------
atom 1 = Any[0.999898, 0.999898, 0.999898]
check lower bound  = -0.04061405485307468
check inequality 1 = 0.9998976341583888
check inequality 2 = 0.9998976341482999
check inequality 3 = 0.999897634147857
check inequality 4 = 0.0003469021566879782
check inequality 5 = 2.9993858363453922
check equality 1 = -0.0003070661102261285
####################################
Solution = Any[0.999898, 0.999898, 0.999898]
####################################
------------------------------------
L0 = 3.039999891198467
====================================
 18.072337 seconds (8.87 M allocations: 436.861 MiB, 2.89% gc time)
3-element Array{Any,1}:
 0.9998976341583888
 0.9998976341482999
 0.999897634147857 
 """