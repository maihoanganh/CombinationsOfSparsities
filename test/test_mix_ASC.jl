include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    eps=1e-2
    k=2 # relaxed order
    r=4 #sparse order

    @polyvar x1 x2 x3 # variables
    x=[x1;x2;x3]

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

    sol = CombinationsOfSparsities.mix_sparsity_adding_spherical_constraints(x,g,h,I,J,W,eps,k,r)
end


"""
Determine L0:
------------------------------------
WARNING: replacing module CombinationsOfSparsities.

The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0006, 0x000b, 0x000d, 0x000f, 0x0016, 0x0018, 0x001a, 0x001c], [0x0002, 0x0007, 0x0009, 0x0010, 0x0012, 0x0014], [0x0003, 0x0008, 0x000a, 0x0011, 0x0013, 0x0015], [0x0005, 0x000c, 0x000e, 0x0017, 0x0019, 0x001b]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001, 0x0004, 0x0006, 0x000b, 0x000d, 0x000f], [0x0002, 0x0007, 0x0009, 0x0010, 0x0012, 0x0014], [0x0003, 0x0008, 0x000a, 0x0011, 0x0013, 0x0015], [0x0005, 0x000c, 0x000e]]]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[]
=============================
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0005, 0x0006, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f, 0x0016, 0x0017, 0x0018, 0x0019, 0x001a, 0x001b, 0x001c], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a, 0x0010, 0x0011, 0x0012, 0x0013, 0x0014, 0x0015]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[]
-----------------------------
block_psi = Array{Array{UInt16,1},1}[[[0x0001, 0x0004, 0x0005, 0x0006, 0x000b, 0x000c, 0x000d, 0x000e, 0x000f], [0x0002, 0x0003, 0x0007, 0x0008, 0x0009, 0x000a, 0x0010, 0x0011, 0x0012, 0x0013, 0x0014, 0x0015]]]
=============================
Termination status = OPTIMAL
Optimal value = 1.062393848838014
====================================
Atom(I[1]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [-3.950451160908658e-17, 0.7021224139930896]
check lower bound  = 0
check inequality 1 = 0.5070241157685165
####################################
Solution = [-3.950451160908658e-17, 0.7021224139930896]
####################################
------------------------------------
atom 2 = [-4.7381590518839786e-17, -0.7021224139930897]
check lower bound  = 0
check inequality 1 = 0.5070241157685162
####################################
Solution = [-4.7381590518839786e-17, -0.7021224139930897]
####################################
====================================
Atom(I[2]):
====================================
Dimension of the null space of Gram matrix = 4
------------------------------------
L0 = 1.062393848838014
====================================
Array{Any,2}
Determine omega0[1]:
------------------------------------
The block-closure operation is stable at sparse order 1!
block_sigma0 = Array{UInt16,1}[[0x0001, 0x0004, 0x0006], [0x0002], [0x0003], [0x0005]]
-----------------------------
block_sigma = Array{Array{UInt16,1},1}[[[0x0001], [0x0002], [0x0003]], [[0x0001], [0x0002], [0x0003]]]
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
Termination status = OPTIMAL
Optimal value = 0.3518364647207896
====================================
Atom(I[1]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [0.0, -0.5931759884037682]
check lower bound  = 0
check inequality 1 = 0.6481422467812127
check inequality 2 = 0.7105360956192267
####################################
Solution = [0.0, -0.5931759884037682]
####################################
------------------------------------
atom 2 = [0.0, 0.5931759884037682]
check lower bound  = 0
check inequality 1 = 0.6481422467812127
check inequality 2 = 0.7105360956192267
####################################
Solution = [0.0, 0.5931759884037682]
####################################
====================================
Atom(I[2]):
====================================
Dimension of the null space of Gram matrix = 2
------------------------------------
atom 1 = [-0.5931957870142465, -0.8429798917872675]
check lower bound  = 0
check inequality 1 = -0.00010249085111047052
check equality 1 = 5.2120345932493706e-5
####################################
Solution = [-0.5931957870142465, -0.8429798917872675]
####################################
------------------------------------
atom 2 = [0.5931957870142465, 0.8429798917872636]
check lower bound  = 0
check inequality 1 = -0.00010249085110380918
check equality 1 = 5.212034593016224e-5
####################################
Solution = [0.5931957870142465, 0.8429798917872636]
####################################
------------------------------------
omega0[1] = 0.35183647
====================================
 25.566870 seconds (9.83 M allocations: 4.751 GiB, 4.45% gc time)
2-element Array{Array{Float64,1},1}:
 [0.0, 0.5931759884037682]               
 [0.5931957870142465, 0.8429798917872636]   
"""