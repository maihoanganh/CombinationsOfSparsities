include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
    eps=1e-2
    k=2 # relaxed order
    r=2 # sparse order

    @polyvar x1 x2 x3 x4# define polnomial variables
    x=[x1;x2;x3;x4] 

    g = [x1^2+x2^2+x3^2+x4^2-2*x1*x3-2*x2*x4-1] #inequality constraints
    h = [x1^2+x2^2-1;x3^2+x4^2-1] #equality constraints

    sol=CombinationsOfSparsities.term_sparsity_adding_spherical_constraints(x,g,h,eps,k,r)
end

"""
Determine L0:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 6, 8, 9, 11, 13, 15], [2, 4, 16, 18, 20, 22, 23, 25, 27, 30, 32, 34], [3, 5, 17, 19, 21, 24, 26, 28, 29, 31, 33, 35], [7, 10, 12, 14]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0006, 0x0008, 0x0009, 0x000b, 0x000d, 0x000f], [0x0002, 0x0004], [0x0003, 0x0005], [0x0007, 0x000a, 0x000c, 0x000e]], Array{UInt16,1}[[0x0001, 0x0006, 0x0008, 0x0009, 0x000b, 0x000d, 0x000f], [0x0002, 0x0004], [0x0003, 0x0005], [0x0007, 0x000a, 0x000c, 0x000e]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0006, 0x0008, 0x0009, 0x000b, 0x000d, 0x000f], [0x0002, 0x0004], [0x0003, 0x0005], [0x0007, 0x000a, 0x000c, 0x000e]], Array{UInt16,1}[[0x0001, 0x0006, 0x0008, 0x0009, 0x000b, 0x000d, 0x000f], [0x0002, 0x0004], [0x0003, 0x0005], [0x0007, 0x000a, 0x000c, 0x000e]]]
Termination status = OPTIMAL
Optimal value = 2.0299999895135175
Rank of moment matrix = 23
------------------------------------
L0 = 2.0299999895135175
====================================
Determine omega0:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 6, 8, 9, 11, 13, 15], [2, 4], [3, 5], [7, 10, 12, 14]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]], Array{UInt16,1}[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]], Array{UInt16,1}[[0x0001], [0x0002, 0x0004], [0x0003, 0x0005]]]
Termination status = OPTIMAL
Optimal value = 1.9999999980162433
Rank of moment matrix = 13
------------------------------------
omega0 = 1.9999999980162433
====================================
Determine omega1:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 2, 4, 6, 8, 9, 11, 13, 15], [3, 5, 7, 10, 12, 14]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0004], [0x0003, 0x0005]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0004], [0x0003, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0004], [0x0003, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0004], [0x0003, 0x0005]]]
Termination status = OPTIMAL
Optimal value = 1.0000000040787809
Rank of moment matrix = 5
------------------------------------
omega1 = 1.0000000040787809
====================================
Determine omega2:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
Termination status = SLOW_PROGRESS
Optimal value = 2.9992794445308477
Rank of moment matrix = 5
------------------------------------
omega2 = 2.9992794445308477
====================================
Determine omega3:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
Termination status = SLOW_PROGRESS
Optimal value = 2.0013557417036205
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[1.0, 0.000358487, 0.500876, 0.299668]
check lower bound  = -0.6624252306665657
check inequality 1 = -0.6612865605509985
check equality 1 = 4.555133629624564e-6
check equality 2 = -0.6593221548589443
check equality 3 = 0.659317597741558
check equality 4 = 0.6593220304199316
check equality 5 = 0.659314018160188
------------------------------------
atom 2 = Any[1.00007, 0.00551439, -0.48088, 0.297152]
check lower bound  = 1.2801258076790512
check inequality 1 = 1.2782749287333592
check equality 1 = 0.00017722776148909603
check equality 2 = -0.680455488227812
check equality 3 = 0.6802782584825662
check equality 4 = 0.6804250783897487
check equality 5 = 0.6905864906543435
------------------------------------
omega3 = 2.0013557417036205
====================================
Determine omega4:
------------------------------------
block_sigma0 = Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]]
-----------------------------
block_sigma = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
-----------------------------
block_psi = Any[Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]], Array{UInt16,1}[[0x0001, 0x0002, 0x0003, 0x0004, 0x0005]]]
Termination status = SLOW_PROGRESS
Optimal value = 4.069272866340988
Rank of moment matrix = 2
------------------------------------
atom 1 = Any[1.00005, 0.000343726, 0.499367, -0.866398]
check lower bound  = 0.66363619134997
check inequality 1 = 0.0019247202413530395
check equality 1 = 0.00010116712262764338
check equality 2 = 1.2437186863278882e-5
check equality 3 = -0.00011360629324741467
check equality 4 = -1.2553808136139821e-5
check equality 5 = -0.0001467080809467447
check equality 6 = -2.3831200126611662e-5
------------------------------------
atom 2 = Any[0.999994, 0.000360638, 0.499312, 0.866388]
check lower bound  = -2.8021199526050165
check inequality 1 = 0.0006853153414605462
check equality 1 = -1.2104037087623531e-5
check equality 2 = -5.9043673409031605e-5
check equality 3 = 7.114572674016273e-5
check equality 4 = 5.8917654844226064e-5
check equality 5 = 7.186868589914575e-5
check equality 6 = 5.162978273642871e-5
------------------------------------
omega4 = 4.069272866340988
====================================
 11.379392 seconds (5.54 M allocations: 284.788 MiB, 3.04% gc time)
0-element Array{Any,1}
"""