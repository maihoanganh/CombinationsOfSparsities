include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
eps=1e-5
k=1 # relaxed order
r=4 #sparse order

@polyvar x1 x2 x3 # variables
x=[x1;x2;x3]

f=Array{Any}(undef, 2)
f[1]=x1*x2
f[2]=x2*x3

f0=f[1]+f[2]
g=[1-x1^2-x2^2]
h=[x2*x3-1/2]



I=Array{Any}(undef, 2)
I[1]=[1;2]
I[2]=[2;3]

J=Array{Any}(undef, 2)
J[1]=[1]
J[2]=[]

W=Array{Any}(undef, 2)
W[1]=[]
W[2]=[1]

opt_val, sol = CombinationsOfSparsities.PuVa_mix_hierarchy(x,f,f0,g,I,J,W,h,eps,k,r)
end


"""
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 4, 5, 6, 11, 12, 13, 14, 15], [2, 3, 7, 8, 9, 10]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]]
-----------------------------
block_psi = Any[]
=============================
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{Int64,1}[[1, 4, 5, 6, 11, 12, 13, 14, 15], [2, 3, 7, 8, 9, 10]]
-----------------------------
block_sigma = Any[]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]]
=============================
WARNING: replacing module CombinationsOfSparsities.

Termination status = OPTIMAL
Optimal value = 7.99914990288441e-5
====================================
Atom(I[1]):
====================================
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[0.7070960264298758, -0.7071177438277267]
check lower bound  = 0
check inequality 1 = -2.942289341056892e-7
####################################
Solution = Any[0.7070960264298758, -0.7071177438277267]
####################################
------------------------------------
atom 2 = Any[-0.707096026429876, 0.7071177438277269]
check lower bound  = 0
check inequality 1 = -2.942289345497784e-7
####################################
Solution = Any[-0.707096026429876, 0.7071177438277269]
####################################
====================================
Atom(I[2]):
====================================
Rank of moment matrix = 4
------------------------------------
atom 1 = Any[0.04666385995922851, 10.203339893666406]
check lower bound  = 0
check equality 1 = -0.02387277608554128
------------------------------------
atom 2 = Any[0.7058346476324053, 0.7088753005548294]
check lower bound  = 0
check equality 1 = 0.0003487479824334505
####################################
Solution = Any[0.7058346476324053, 0.7088753005548294]
####################################
------------------------------------
atom 3 = Any[-0.7064967133755139, -0.7079362254585579]
check lower bound  = 0
check equality 1 = 0.00015461656593795947
####################################
Solution = Any[-0.7064967133755139, -0.7079362254585579]
####################################
------------------------------------
atom 4 = Any[-0.04600179421073918, -10.204278968765447]
check lower bound  = 0
check equality 1 = -0.0305848588098781
 25.655892 seconds (3.64 M allocations: 216.662 MiB, 13.23% gc time)
(7.99914990288441e-5, Any[[0.7070960264298758 -0.707096026429876; -0.7071177438277267 0.7071177438277269], [0.7058346476324053 -0.7064967133755139; 0.7088753005548294 -0.7079362254585579]])
"""