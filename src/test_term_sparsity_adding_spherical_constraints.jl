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