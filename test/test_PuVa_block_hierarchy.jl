include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    eps=1e-7
    k=2 # relaxed order
    r=10 # sparse order

    @polyvar x1 x2 x3 # define polnomial variables
    x=[x1;x2;x3]

    f = x1+x2+x3 # objective function
    g = [x1;x2;x3] # inequality constraints
    h = [x1*x2*x3-1] # equality constraints

    opt_val,sol = CombinationsOfSparsities.PuVa_block_hierarchy(x,f,g,h,eps,k,r,knownlb=0)

    g=[g;opt_val-f]

    eps=1e-7
    k=1 # relaxed order
    r=2 # sparse order

    sol=CombinationsOfSparsities.block_ASC(x,g,h,eps,k,r) # compute minimizers
end


"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]], [[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]], [[0x0001, 0x0002, 0x0003, 0x0004, 0x0005, 0x0006, 0x0007, 0x0008, 0x0009, 0x000a]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001, 0x0002, 0x0003, 0x0004]]]
WARNING: replacing module CombinationsOfSparsities.
Termination status = OPTIMAL
Optimal value = 3.0000002895711915
Dimension of the null space of Gram matrix = 11
Determine L0:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001, 0x0002, 0x0003, 0x0004]], [[0x0001, 0x0002, 0x0003, 0x0004]], [[0x0001, 0x0002, 0x0003, 0x0004]], [[0x0001, 0x0002, 0x0003, 0x0004]], [[0x0001, 0x0002, 0x0003, 0x0004]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]]]
Termination status = OPTIMAL
Optimal value = 3.0000003650799525
Dimension of the null space of Gram matrix = 1
------------------------------------
atom 1 = [1.0000203491343966, 1.0000203492526776, 1.0000203499704754]
check lower bound  = 0.00012173287744721151
check inequality 1 = 1.0000203491343966
check inequality 2 = 1.0000203492526776
check inequality 3 = 1.0000203499704754
check inequality 4 = -6.075878635813936e-5
check inequality 5 = 3.0001220979573997
check equality 1 = 6.104959985875347e-5
####################################
Solution = [1.0000203491343966, 1.0000203492526776, 1.0000203499704754]
####################################
------------------------------------
L0 = 3.0000003650799525
====================================
 17.472473 seconds (6.97 M allocations: 335.934 MiB, 1.75% gc time)
1-element Array{Array{Float64,1},1}:
 [1.0000203491343966, 1.0000203492526776, 1.0000203499704754]
 """