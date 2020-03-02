include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials
@time begin
    k=2 # relaxed order
    r=2 # sparse order

    @polyvar x1 x2 # define polnomial variables
    x=[x1;x2]

    f = x1^4+x2^4-x1*x2+1 # objective function
    g = [1-x1^2-2*x2^2;x1^2+x2^2-1] # inequality constraints
    h = [] # equality constraints

    opt_val, sol = CombinationsOfSparsities.Pu_block_hierarchy(x,f,g,h,k,r)
end


"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 4, 5, 6], [2, 3]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0003]], Array{UInt16,1}[[0x0001], [0x0002, 0x0003]]]
-----------------------------
block_psi = Any[]
Termination status = SLOW_PROGRESS
Optimal value = 1.9999885838473876
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[0.9999999872200782, 1.8653288232151585e-5]
check lower bound  = -7.28825506768338e-6
check inequality 1 = 2.486395311862566e-8
check inequality 2 = -2.5211898346810813e-8
####################################
Solution = Any[0.9999999872200782, 1.8653288232151585e-5]
####################################
------------------------------------
atom 2 = Any[-0.9999999872200781, -1.8653288232151585e-5]
check lower bound  = -7.288255068127469e-6
check inequality 1 = 2.4863953340670264e-8
check inequality 2 = -2.5211898568855418e-8
####################################
Solution = Any[-0.9999999872200781, -1.8653288232151585e-5]
####################################
  4.933601 seconds (2.59 M allocations: 131.997 MiB, 2.07% gc time)
(1.9999885838473876, [0.9999999872200782 -0.9999999872200781; 1.8653288232151585e-5 -1.8653288232151585e-5])
"""