include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    eps=1e-7
    k=1 # relaxed order
    r=10 # sparse order

    @polyvar x1 x2 x3 x4# define polnomial variables
    x=[x1;x2;x3;x4] 

    g = [x1^2+x2^2+x3^2+x4^2-2*x1*x3-2*x2*x4-1.0] #inequality constraints
    h = [x1^2+x2^2-1;x3^2+x4^2-1] #equality constraints

    sol=CombinationsOfSparsities.block_ASC(x,g,h,eps,k,r)
end

"""
Determine L0:
------------------------------------
The block-closure operation is stable at sparse order 2!
block_sigma0 = Array{Int64,1}[[1, 6, 8, 9, 11, 13, 15], [2, 4], [3, 5], [7, 10, 12, 14]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]], [[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]
WARNING: replacing module CombinationsOfSparsities.
]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]], [[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]]]
Termination status = OPTIMAL
Optimal value = 2.0000002926768916
Dimension of the null space of Gram matrix = 13
------------------------------------
L0 = 2.0000002926768916
====================================
Determine omega0:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1], [2, 4], [3, 5]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]]]
Termination status = OPTIMAL
Optimal value = 1.999999956856609
Dimension of the null space of Gram matrix = 5
------------------------------------
omega0 = 1.999999956856609
====================================
Determine omega1:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 4], [3, 5]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]], [[0x0001]]]
Termination status = SLOW_PROGRESS
Optimal value = 1.0000001068472149
Dimension of the null space of Gram matrix = 3
------------------------------------
omega1 = 1.0000001068472149
====================================
Determine omega2:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]]]
Termination status = SLOW_PROGRESS
Optimal value = 2.999020891000328
Dimension of the null space of Gram matrix = 1
------------------------------------
atom 1 = [1.0000008021563165, 0.0006137650000046859, -0.3252888930722904, -0.00019965463590494666]
check lower bound  = -0.894433536161416
check inequality 1 = 0.7563934379306378
check equality 1 = 1.981020751662399e-6
check equality 2 = -0.8941870961818303
check equality 3 = 0.8941850720176876
check equality 4 = 0.8941868263209265
------------------------------------
omega2 = 2.999020891000328
====================================
Determine omega3:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]]]
Termination status = SLOW_PROGRESS
Optimal value = 1.999402753207286
Dimension of the null space of Gram matrix = 1
------------------------------------
atom 1 = [0.9999962692194401, 0.0005886122850475202, 0.495203643250304, -0.855044773307536]
check lower bound  = -0.013488942141763216
check inequality 1 = -0.013075914226557295
check equality 1 = -7.115082778952342e-6
check equality 2 = -0.02367178735108999
check equality 3 = 0.02367885929047797
check equality 4 = 0.023671547719964048
check equality 5 = 0.023877018004291983
------------------------------------
omega3 = 1.999402753207286
====================================
Determine omega4:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]], [[0x0001]]]
Termination status = SLOW_PROGRESS
Optimal value = 4.093301691912711
Dimension of the null space of Gram matrix = 0
------------------------------------
omega4 = 4.093301691912711
====================================
 15.479054 seconds (5.97 M allocations: 285.470 MiB, 1.59% gc time)
0-element Array{Array{Float64,1},1}
"""