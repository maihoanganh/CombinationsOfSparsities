include("./CombinationsOfSparsities.jl")
using .CombinationsOfSparsities

k=2 # relaxed order
r=1 # sparse order

@polyvar x1 x2 # define polnomial variables
x=[x1;x2]

f = x1^4+x2^4-x1*x2 # objective function
g = [x1;1-x1^2-x2^2] # inequality constraints
h = [x1+1/2-(x2+1/2)^3] # equality constraints

opt_val, sol = Pu_block_hierarchy(x,f,g,h,k,r)
