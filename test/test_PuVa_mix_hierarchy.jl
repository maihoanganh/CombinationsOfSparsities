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
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[]
=============================
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001, 0x0004, 0x0005, 0x0006], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a]]]
=============================
WARNING: replacing module CombinationsOfSparsities.
Termination status = OPTIMAL
Optimal value = 7.99914990288441e-5
====================================
Atom(I[1]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [0.7070924156516603, -0.7071141329386108]
check lower bound  = 0
check inequality 1 = 9.91872637634561e-6
####################################
Solution = [0.7070924156516603, -0.7071141329386108]
####################################
------------------------------------
atom 2 = [-0.7070924156516599, 0.7071141329386114]
check lower bound  = 0
check inequality 1 = 9.918726376012543e-6
####################################
Solution = [-0.7070924156516599, 0.7071141329386114]
####################################
====================================
Atom(I[2]):
====================================
Dimension of the null space of Gram matrix = 4
------------------------------------
atom 1 = [-0.04659137774972165, -10.214970171130298]
check lower bound  = 0
check equality 1 = -0.024070466054729456
------------------------------------
atom 2 = [-0.7059163897012568, -0.7081746164819651]
check lower bound  = 0
check equality 1 = -8.793145497903998e-5
####################################
Solution = [-0.7059163897012568, -0.7081746164819651]
####################################
------------------------------------
atom 3 = [0.7064550437866567, 0.70767550575466]
check lower bound  = 0
check equality 1 = -5.906959534723866e-5
####################################
Solution = [0.7064550437866567, 0.70767550575466]
####################################
------------------------------------
atom 4 = [0.046052723666656865, 10.215469281886337]
check lower bound  = 0
check equality 1 = -0.029549816036066878
  9.898107 seconds (3.82 M allocations: 295.371 MiB, 2.33% gc time)
(7.99914990288441e-5, Array{Float64,1}[[-0.7070924156516599, 0.7071141329386114], [0.7064550437866567, 0.70767550575466]])
"""