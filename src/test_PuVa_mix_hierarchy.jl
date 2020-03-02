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
