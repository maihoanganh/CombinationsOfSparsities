function Pu_sparse_hierarchy(x,f0,g,h,I,J,W,k)

    n0=length(x) # Number of variables
    m0=length(g) # Number of constraints
    l0=length(h) # Number of constraints

    supp_f0,coe_f0=info(f0,x,n0)
    lf0=length(coe_f0)

    supp_g=Array{Any}(undef, m0)
    coe_g=Array{Any}(undef, m0)
    lg=Array{Any}(undef, m0)
    for j=1:m0
        supp_g[j],coe_g[j]=info(g[j],x,n0)
        lg[j]=length(coe_g[j])
    end

    supp_h=Array{Any}(undef, l0)
    coe_h=Array{Any}(undef, l0)
    lh=Array{Any}(undef, l0)
    for j=1:l0
        supp_h[j],coe_h[j]=info(h[j],x,n0)
        lh[j]=length(coe_h[j])
    end

    p=length(I)

    n=Array{Int64}(undef, p)
    m=Array{UInt16}(undef, p)
    l=Array{UInt16}(undef, p)

    u=Array{Any}(undef, p)
    w=Array{Any}(undef, p)

    lu0=Array{Any}(undef, p)
    lu=Array{Any}(undef, p)
    lw=Array{Any}(undef, p)

    basis_sigma0=Array{Any}(undef, p)
    basis_sigma=Array{Any}(undef, p)
    basis_psi=Array{Any}(undef, p)

    largest_supp=[]

    #Degrees of constraints
    for t=1:p
        n[t]=length(I[t])

        lu0[t]=binomial(n[t]+k,n[t])
        basis_sigma0[t]=zeros(n0,lu0[t])
        basis_sigma0[t][I[t],:]=get_basis(n[t],k)

        for i=1:lu0[t]
            for r=i:lu0[t]
                nota=basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                if largest_supp==[]
                    largest_supp=nota
                else
                    largest_supp=[largest_supp nota]
                end
            end
        end

        m[t]=length(J[t])
        u[t]=Array{Int16}(undef, m[t])
        lu[t]=Array{UInt16}(undef, m[t])

        basis_sigma[t]=Array{Any}(undef, m[t])
        for j=1:m[t]
            u[t][j]=ceil(Int16,.5*maxdegree(g[J[t][j]]))
            lu[t][j]=binomial(n[t]+k-u[t][j],n[t])
            basis_sigma[t][j]=basis_sigma0[t][:,1:lu[t][j]]

            for i=1:lu[t][j]
                for r=i:lu[t][j]
                    for z=1:lg[J[t][j]]
                        nota=supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                        largest_supp=[largest_supp nota]
                    end
                end
            end
        end

        l[t]=length(W[t])
        w[t]=Array{Int16}(undef, l[t])
        lw[t]=Array{UInt16}(undef, l[t])

        basis_psi[t]=Array{Any}(undef, l[t])
        for j=1:l[t]
            w[t][j]=ceil(Int16,.5*maxdegree(h[W[t][j]]))
            lw[t][j]=binomial(n[t]+k-w[t][j],n[t])
            basis_psi[t][j]=basis_sigma0[t][:,1:lw[t][j]]

            for i=1:lw[t][j]
                for r=i:lw[t][j]
                    for z=1:lh[W[t][j]]
                        nota=supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                        largest_supp=[largest_supp nota]
                    end
                end
            end
        end

    end


    largest_supp=unique(largest_supp,dims=2)
    largest_supp=sortslices(largest_supp,dims=2)
    llargest_supp=size(largest_supp,2)

    model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=[AffExpr(0) for i=1:llargest_supp]

    G0=Array{Any}(undef, p)
    G=Array{Any}(undef, p)
    H=Array{Any}(undef, p)

    for t=1:p
        G0[t]=@variable(model, [1:lu0[t], 1:lu0[t]],PSD)
        for i=1:lu0[t]
            for r=i:lu0[t]
                nota=basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                if i==r
                  cons[Locb]+=G0[t][i,r]
                else
                  cons[Locb]+=2*G0[t][i,r]
                end
            end
        end

        G[t]=Array{Any}(undef, m[t])
        for j=1:m[t]
            G[t][j]=@variable(model, [1:lu[t][j], 1:lu[t][j]],PSD)
            for i=1:lu[t][j]
                for r=i:lu[t][j]
                    for z=1:lg[J[t][j]]
                        nota=supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        if i==r
                          cons[Locb]+=coe_g[J[t][j]][z]*G[t][j][i,r]
                        else
                          cons[Locb]+=2*coe_g[J[t][j]][z]*G[t][j][i,r]
                        end
                    end
                end
            end
        end

        H[t]=Array{Any}(undef, l[t])
        for j=1:l[t]
            H[t][j]=@variable(model, [1:lw[t][j], 1:lw[t][j]],Symmetric)
            for i=1:lw[t][j]
                for r=i:lw[t][j]
                    for z=1:lh[W[t][j]]
                        nota=supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        if i==r
                          cons[Locb]+=coe_h[W[t][j]][z]*H[t][j][i,r]
                        else
                          cons[Locb]+=2*coe_h[W[t][j]][z]*H[t][j][i,r]
                        end
                    end
                end
            end
        end
    end

    bc=zeros(1,llargest_supp)
    for i=1:lf0
        Locb=bfind(largest_supp,llargest_supp,supp_f0[:,i],n0)
        bc[Locb]=coe_f0[i]
    end
    @constraint(model, cons[2:end].==bc[2:end])
    @variable(model, lambda)
    @constraint(model, cons[1]+lambda==bc[1])
    @objective(model, Max, lambda)
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)

    #extraction of optimizers
    sol=Array{Any}(undef,2)
    for t=1:p
        println("====================================")
        println("Atom(I[",t,"]):")
        println("====================================")
        sol[t]=extract_optimizers(value.(G0[t]),lu0[t],basis_sigma0[t],n0,n[t],m[t],l[t],opt_val,0,g[J[t]],h[W[t]],x,I[t])
    end

    return opt_val,sol

end
