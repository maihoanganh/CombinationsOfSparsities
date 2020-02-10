function Pu_mix_hierarchy(x,f,f0,g,I,J,W,h,k,r)
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

    supp_f=Array{Any}(undef, p)
    coe_f=Array{Any}(undef, p)
    lf=Array{Any}(undef, p)

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

    supp_U=Array{Any}(undef, p)
    lsupp_U=Array{Any}(undef, p)

    block_sigma0=Array{Any}(undef, p)
    block_sigma=Array{Any}(undef, p)
    block_psi=Array{Any}(undef, p)

    lblock_sigma0=Array{Any}(undef, p)
    lblock_sigma=Array{Any}(undef, p)
    lblock_psi=Array{Any}(undef, p)

    lt_block_sigma0=Array{Any}(undef, p)
    lt_block_sigma=Array{Any}(undef, p)
    lt_block_psi=Array{Any}(undef, p)

    largest_supp=[]
    #Degrees of constraints
    for t=1:p

        n[t]=length(I[t])

        supp_f[t],coe_f[t]=info(f[t],x,n0)
        lf[t]=length(coe_f[t])

        A=supp_f[t]

        lu0[t]=binomial(n[t]+k,n[t])
        basis_sigma0[t]=zeros(n0,lu0[t])
        basis_sigma0[t][I[t],:]=get_basis(n[t],k)


        m[t]=length(J[t])
        u[t]=Array{Int16}(undef, m[t])
        lu[t]=Array{UInt16}(undef, m[t])

        basis_sigma[t]=Array{Any}(undef, m[t])
        for j=1:m[t]
            u[t][j]=ceil(Int16,.5*maxdegree(g[J[t][j]]))
            lu[t][j]=binomial(n[t]+k-u[t][j],n[t])
            basis_sigma[t][j]=basis_sigma0[t][:,1:lu[t][j]]

            A=[A supp_g[J[t][j]]]
        end

        l[t]=length(W[t])
        w[t]=Array{Int16}(undef, l[t])
        lw[t]=Array{UInt16}(undef, l[t])

        basis_psi[t]=Array{Any}(undef, l[t])
        for j=1:l[t]
            w[t][j]=ceil(Int16,.5*maxdegree(h[W[t][j]]))
            lw[t][j]=binomial(n[t]+k-w[t][j],n[t])
            basis_psi[t][j]=basis_sigma0[t][:,1:lw[t][j]]

            A=[A supp_h[W[t][j]]]
        end

        supp_U[t]=[A zeros(n0,1)]

        supp_U[t],lsupp_U[t],block_sigma0[t],block_sigma[t],block_psi[t],lblock_sigma0[t],lblock_sigma[t],lblock_psi[t],lt_block_sigma0[t],lt_block_sigma[t],lt_block_psi[t]=get_blocks(r,n0,m[t],l[t],supp_U[t],lu0[t],lu[t],lw[t],lg[J[t]],lh[W[t]],supp_g[J[t]],supp_h[W[t]],coe_g[J[t]],coe_h[W[t]],basis_sigma0[t],basis_sigma[t],basis_psi[t])

        println("block_sigma0 = ",block_sigma0[t])
        println("-----------------------------")
        println("block_sigma = ",block_sigma[t])
        println("-----------------------------")
        println("block_psi = ",block_psi[t])
        println("=============================")
        if t==1
            largest_supp=supp_U[t]
        else
            largest_supp=[largest_supp supp_U[t]]
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
        G0[t]=Array{Any}(undef, lblock_sigma0[t])
        G[t]=Array{Any}(undef, m[t])
        H[t]=Array{Any}(undef, l[t])

        for j=1:lblock_sigma0[t]
            if lt_block_sigma0[t][j]==1
                G0[t][j]=@variable(model, lower_bound=0)
                nota=UInt8(2)*basis_sigma0[t][:,block_sigma0[t][j]]
                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                cons[Locb]+=G0[t][j]
            else
                G0[t][j]=@variable(model, [1:lt_block_sigma0[t][j], 1:lt_block_sigma0[t][j]],PSD)
                for p=1:lt_block_sigma0[t][j]
                    for q=p:lt_block_sigma0[t][j]
                        nota=basis_sigma0[t][:,block_sigma0[t][j][p]]+basis_sigma0[t][:,block_sigma0[t][j][q]]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        if p==q
                            cons[Locb]+=G0[t][j][p,q]
                        else
                            cons[Locb]+=2*G0[t][j][p,q]
                        end
                    end
                end
            end
        end


        for i=1:m[t]
            G[t][i]=Array{Any}(undef, lblock_sigma[t][i])
            for j=1:lblock_sigma[t][i]
                if lt_block_sigma[t][i][j]==1
                    G[t][i][j]=@variable(model, lower_bound=0)
                    for z=1:lg[J[t][i]]
                        nota=supp_g[J[t][i]][:,z]+UInt8(2)*basis_sigma[t][i][:,block_sigma[t][i][j]]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        cons[Locb]+=coe_g[J[t][i]][z]*G[t][i][j]
                    end
                else
                    G[t][i][j]=@variable(model, [1:lt_block_sigma[t][i][j], 1:lt_block_sigma[t][i][j]],PSD)
                    for p=1:lt_block_sigma[t][i][j]
                        for q=p:lt_block_sigma[t][i][j]
                            for z=1:lg[J[t][i]]
                                nota=basis_sigma[t][i][:,block_sigma[t][i][j][p]]+basis_sigma[t][i][:,block_sigma[t][i][j][q]]+supp_g[J[t][i]][:,z]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if p==q
                                  cons[Locb]+=coe_g[J[t][i]][z]*G[t][i][j][p,q]
                                else
                                  cons[Locb]+=2*coe_g[J[t][i]][z]*G[t][i][j][p,q]
                                end
                            end
                        end
                    end
                end
            end
        end

        for i=1:l[t]
            H[t][i]=Array{Any}(undef, lblock_psi[t][i])
            for j=1:lblock_psi[t][i]
                if lt_block_psi[t][i][j]==1
                    H[t][i][j]=@variable(model)
                    for z=1:lh[W[t][i]]
                        nota=supp_h[W[t][i]][:,z]+UInt8(2)*basis_psi[t][i][:,block_psi[t][i][j]]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        cons[Locb]+=coe_h[W[t][i]][z]*H[t][i][j]
                    end
                else
                    H[t][i][j]=@variable(model, [1:lt_block_psi[t][i][j], 1:lt_block_psi[t][i][j]],Symmetric)
                    for p=1:lt_block_psi[t][i][j]
                        for q=p:lt_block_psi[t][i][j]
                            for z=1:lh[W[t][i]]
                                nota=basis_psi[t][i][:,block_psi[t][i][j][p]]+basis_psi[t][i][:,block_psi[t][i][j][q]]+supp_h[W[t][i]][:,z]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if p==q
                                  cons[Locb]+=coe_h[W[t][i]][z]*H[t][i][j][p,q]
                                else
                                  cons[Locb]+=2*coe_h[W[t][i]][z]*H[t][i][j][p,q]
                                end
                            end
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
        Gr=zeros(lu0[t],lu0[t])
        for j=1:lblock_sigma0[t]
            if lt_block_sigma0[t][j]==1
                Gr[block_sigma0[t][j][1],block_sigma0[t][j][1]]=value(G0[t][j])
            else
                Gr[block_sigma0[t][j],block_sigma0[t][j]]=value.(G0[t][j])
            end
        end
        sol[t]=extract_optimizers(Gr,lu0[t],basis_sigma0[t],n0,n[t],m[t],l[t],opt_val,0,g[J[t]],h[W[t]],x,I[t])
    end

    return opt_val, sol

end
