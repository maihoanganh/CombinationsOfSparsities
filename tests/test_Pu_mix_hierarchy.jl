include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

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

"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 4, 5, 6], [2, 3]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0003]]]
-----------------------------
block_psi = Any[]
=============================
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 4, 5, 6], [2, 3]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0003]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0003]]]
=============================
Termination status = OPTIMAL
Optimal value = 5.012232764705938e-9
====================================
Atom(I[1]):
====================================
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[0.707106, -0.70711]
check lower bound  = 0
check inequality 1 = -4.031208835630906e-6
####################################
Solution = Any[0.707106, -0.70711]
####################################
------------------------------------
atom 2 = Any[-0.707106, 0.70711]
check lower bound  = 0
check inequality 1 = -4.031208835630906e-6
####################################
Solution = Any[-0.707106, 0.70711]
####################################
====================================
Atom(I[2]):
====================================
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[0.707106, 0.707124]
check lower bound  = 0
check inequality 1 = -2.2885817929418906e-5
check equality 1 = 1.144274092490516e-5
####################################
Solution = Any[0.707106, 0.707124]
####################################
------------------------------------
atom 2 = Any[-0.707106, -0.707124]
check lower bound  = 0
check inequality 1 = -2.2885817929418906e-5
check equality 1 = 1.144274092490516e-5
####################################
Solution = Any[-0.707106, -0.707124]
####################################
  8.603841 seconds (2.72 M allocations: 140.270 MiB, 13.91% gc time)
(5.012232764705938e-9, Any[[0.707106 -0.707106; -0.70711 0.70711], [0.707106 -0.707106; 0.707124 -0.707124]])
"""