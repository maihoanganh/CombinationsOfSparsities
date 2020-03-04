include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    k=2 # relaxed order
    r=10 #sparse order

    @polyvar x1 x2 x3 # variables
    x=[x1;x2;x3]

    f=Array{Any}(undef, 2)
    f[1]=x1*x2
    f[2]=x2*x3

    f0=f[1]+f[2]
    g=[1-x1^2-x2^2;
       1-x2^2-x3^2]
    h=[x2*x3-1/2]



    I=Array{Any}(undef, 2)
    I[1]=[1;2]
    I[2]=[2;3]

    J=Array{Any}(undef, 2)
    J[1]=[1]
    J[2]=[2]

    W=Array{Any}(undef, 2)
    W[1]=[]
    W[2]=[1]

    opt_val, sol = CombinationsOfSparsities.Pu_mix_hierarchy(x,f,f0,g,I,J,W,h,k,r)
end

"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001], [0x0002, 0x0003]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[]
=============================
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001], [0x0002, 0x0003]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001], [0x0002, 0x0003]]]
=============================
WARNING: replacing module CombinationsOfSparsities.
Termination status = OPTIMAL
Optimal value = 5.012232764705938e-9
====================================
Atom(I[1]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [-0.7071063667889375, 0.7071099945386363]
check lower bound  = 0
check inequality 1 = -3.958329881648126e-6
####################################
Solution = [-0.7071063667889375, 0.7071099945386363]
####################################
------------------------------------
atom 2 = [0.7071063667889375, -0.7071099945386364]
check lower bound  = 0
check inequality 1 = -3.958329881870171e-6
####################################
Solution = [0.7071063667889375, -0.7071099945386364]
####################################
====================================
Atom(I[2]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [-0.7071028153145639, -0.7071211477140005]
check lower bound  = 0
check inequality 1 = -1.4708970147436773e-5
check equality 1 = 7.354317035357383e-6
####################################
Solution = [-0.7071028153145639, -0.7071211477140005]
####################################
------------------------------------
atom 2 = [0.707102815314564, 0.7071211477140007]
check lower bound  = 0
check inequality 1 = -1.4708970148102907e-5
check equality 1 = 7.354317035579427e-6
####################################
Solution = [0.707102815314564, 0.7071211477140007]
####################################
  4.878905 seconds (2.23 M allocations: 111.555 MiB, 1.54% gc time)
(5.012232764705938e-9, Array{Float64,1}[[0.7071063667889375, -0.7071099945386364], [0.707102815314564, 0.7071211477140007]])
"""