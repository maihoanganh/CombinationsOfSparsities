include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

k=2 # relaxed order

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

opt_val, sol = CombinationsOfSparsities.Pu_sparse_hierarchy(x,f0,g,h,I,J,W,k)

"""
Termination status = OPTIMAL
Optimal value = 8.408748567389877e-9
====================================
Atom(I[1]):
====================================

Rank of moment matrix = 2
------------------------------------
atom 1 = Any[-0.707117, 0.707109]
check lower bound  = 0
check inequality 1 = -1.7941442903790517e-5
####################################
Solution = Any[-0.707117, 0.707109]
####################################
------------------------------------
atom 2 = Any[0.707117, -0.707109]
check lower bound  = 0
check inequality 1 = -1.7941442903568472e-5
####################################
Solution = Any[0.707117, -0.707109]
####################################
====================================
Atom(I[2]):
====================================
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[0.707117, 0.707114]
check lower bound  = 0
check inequality 1 = -2.523046184066402e-5
check equality 1 = 1.2615225078560499e-5
####################################
Solution = Any[0.707117, 0.707114]
####################################
------------------------------------
atom 2 = Any[-0.707117, -0.707114]
check lower bound  = 0
check inequality 1 = -2.523046184021993e-5
check equality 1 = 1.2615225078449477e-5
####################################
Solution = Any[-0.707117, -0.707114]
####################################
 98.158908 seconds (7.76 M allocations: 66.556 GiB, 20.48% gc time)
(8.408748567389877e-9, Any[[-0.707117 0.707117; 0.707109 -0.707109], [0.707117 -0.707117; 0.707114 -0.707114]])
"""