include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

using DynamicPolynomials

@time begin
eps=1e-2
k=2 # relaxed order

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

sol = CombinationsOfSparsities.correlative_sparsity_adding_spherical_constraints(x,g,h,I,J,W,eps,k)
end