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
    have_sol=1
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
        if sol[t]==[]
            have_sol=0
        end
    end
    
    if have_sol==0
        sol=[]
    end
    
    return opt_val, sol

end






function PuVa_mix_hierarchy(x,f,f0,g,I,J,W,h,eps,k,r)
    n0=length(x) # Number of variables
    m0=length(g) # Number of constraints
    l0=length(h) # Number of constraints
    
    d=1+floor(Int16,.5*maxdegree(f0))
    p=length(I)
    
    #Determine hat I_j and s_j
    hat_I=Array{Any}(undef, p)
    s=Array{Any}(undef, p)
    theta=Array{Any}(undef, p)
    hat_theta=Array{Any}(undef, p)
    for j=1:p
        if j>=2
            union_I=0
            for i=1:j-1
                if union_I==0
                    union_I=I[1]
                else
                    union_I=[union_I;I[i]]
                end
            end

            union_I=unique(union_I)
            hat_I[j]=intersect(I[j],union_I)
            hat_theta[j]=sum(x[hat_I[j]].^2)+1
            
            for i=1:j-1
                test_I=intersect(hat_I[j],I[i])
                if length(test_I)==length(hat_I[j])
                    s[j]=i
                    break
                end
            end
        end
        
        theta[j]=sum(x[I[j]].^2)+1
        
    end
    s[2]=1
    
    
    #Define theta_j, hat_theta_j and D_j
    D=Array{Any}(undef, p)
    Theta=Array{Any}(undef, p)
    omega=Array{Any}(undef, p)
    for j=1:p
        if j==1
            D[j]=1
        else
            D[j]=hat_theta[j]
        end
        for i=(j+1):p
            if s[i]==j
                D[j]*=hat_theta[i]
            end
        end
        Theta[j]=D[j]*theta[j]
        omega[j]=ceil(Int16,.5*maxdegree(Theta[j]))
    end
    
    

    psi=0
    phi0=1
    phi=Array{Any}(undef, p)
    
    supp_phik=Array{Any}(undef, p)
    coe_phik=Array{Any}(undef, p)
    lphik=Array{Any}(undef, p)
    for j=1:p
        psi+=theta[j]^d
        phi0 *=Theta[j]

        phi[j]=1
        for r=1:p
            if j!=r
                phi[j]*=Theta[r]
            end
        end
        supp_phik[j],coe_phik[j]=info(phi[j]^k,x,n0)
        lphik[j]=length(coe_phik[j])
    end
    ############
    supp_phik0,coe_phik0=info(phi0^k,x,n0)
    lphik0=length(coe_phik0)    
    
    F=phi0^k*(f0+eps*psi)
    
    supp_F,coe_F=info(F,x,n0)
    lF=length(coe_F)

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

    

    supp_Thetakf=Array{Any}(undef, p)
    coe_Thetakf=Array{Any}(undef, p)
    lThetakf=Array{Any}(undef, p)
    
    supp_Thetak=Array{Any}(undef, p)
    coe_Thetak=Array{Any}(undef, p)
    lThetak=Array{Any}(undef, p)

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

        supp_Thetakf[t],coe_Thetakf[t]=info(Theta[t]^k*(f[t]+eps*theta[t]^d),x,n0)
        lThetakf[t]=length(coe_Thetakf[t])

        A=supp_Thetakf[t]

        lu0[t]=binomial(n[t]+d+k*omega[t],n[t])
        basis_sigma0[t]=zeros(n0,lu0[t])
        basis_sigma0[t][I[t],:]=get_basis(n[t],d+k*omega[t])


        m[t]=length(J[t])
        u[t]=Array{Int16}(undef, m[t])
        lu[t]=Array{UInt16}(undef, m[t])

        basis_sigma[t]=Array{Any}(undef, m[t])
        for j=1:m[t]
            u[t][j]=ceil(Int16,.5*maxdegree(g[J[t][j]]))
            lu[t][j]=binomial(n[t]+d+k*omega[t]-u[t][j],n[t])
            basis_sigma[t][j]=basis_sigma0[t][:,1:lu[t][j]]

            A=[A supp_g[J[t][j]]]
        end

        l[t]=length(W[t])
        w[t]=Array{Int16}(undef, l[t])
        lw[t]=Array{UInt16}(undef, l[t])

        basis_psi[t]=Array{Any}(undef, l[t])
        for j=1:l[t]
            w[t][j]=ceil(Int16,.5*maxdegree(h[W[t][j]]))
            lw[t][j]=binomial(n[t]+d+k*omega[t]-w[t][j],n[t])
            basis_psi[t][j]=basis_sigma0[t][:,1:lw[t][j]]

            A=[A supp_h[W[t][j]]]
        end

        supp_Thetak[t],coe_Thetak[t]=info(Theta[t]^k,x,n0)
        lThetak[t]=length(coe_Thetak[t])
        
        supp_U[t]=[A supp_Thetak[t]]

        supp_U[t],lsupp_U[t],block_sigma0[t],block_sigma[t],block_psi[t],lblock_sigma0[t],lblock_sigma[t],lblock_psi[t],lt_block_sigma0[t],lt_block_sigma[t],lt_block_psi[t]=get_blocks(r,n0,m[t],l[t],supp_U[t],lu0[t],lu[t],lw[t],lg[J[t]],lh[W[t]],supp_g[J[t]],supp_h[W[t]],coe_g[J[t]],coe_h[W[t]],basis_sigma0[t],basis_sigma[t],basis_psi[t])

        println("block_sigma0 = ",block_sigma0[t])
        println("-----------------------------")
        println("block_sigma = ",block_sigma[t])
        println("-----------------------------")
        println("block_psi = ",block_psi[t])
        println("=============================")
    end
    
    
    largest_supp=0
    
    #Degrees of constraints
    for t=1:p
        for j=1:lblock_sigma0[t]
            for p=1:lt_block_sigma0[t][j]
                for q=p:lt_block_sigma0[t][j]
                    for y=1:lphik[t]
                        nota=supp_phik[t][:,y]+basis_sigma0[t][:,block_sigma0[t][j][p]]+basis_sigma0[t][:,block_sigma0[t][j][q]]
                        if largest_supp==0
                            largest_supp=nota
                        else
                            largest_supp=[largest_supp nota]
                        end
                    end
                end
            end
        end


        for i=1:m[t]
            for j=1:lblock_sigma[t][i]
                for p=1:lt_block_sigma[t][i][j]
                    for q=p:lt_block_sigma[t][i][j]
                        for z=1:lg[J[t][i]]
                            for y=1:lphik[t]
                                nota=supp_phik[t][:,y]+basis_sigma[t][i][:,block_sigma[t][i][j][p]]+basis_sigma[t][i][:,block_sigma[t][i][j][q]]+supp_g[J[t][i]][:,z]
                                largest_supp=[largest_supp nota]
                            end
                        end
                    end
                end
            end
        end

        for i=1:l[t]
            for j=1:lblock_psi[t][i]
                for p=1:lt_block_psi[t][i][j]
                    for q=p:lt_block_psi[t][i][j]
                        for z=1:lh[W[t][i]]
                            for y=1:lphik[t]
                                nota=supp_phik[t][:,y]+basis_psi[t][i][:,block_psi[t][i][j][p]]+basis_psi[t][i][:,block_psi[t][i][j][q]]+supp_h[W[t][i]][:,z]
                                largest_supp=[largest_supp nota]
                            end
                        end
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
        G0[t]=Array{Any}(undef, lblock_sigma0[t])
        G[t]=Array{Any}(undef, m[t])
        H[t]=Array{Any}(undef, l[t])

        for j=1:lblock_sigma0[t]
            G0[t][j]=@variable(model, [1:lt_block_sigma0[t][j], 1:lt_block_sigma0[t][j]],PSD)
            for p=1:lt_block_sigma0[t][j]
                for q=p:lt_block_sigma0[t][j]
                    for y=1:lphik[t]
                        nota=supp_phik[t][:,y]+basis_sigma0[t][:,block_sigma0[t][j][p]]+basis_sigma0[t][:,block_sigma0[t][j][q]]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        if p==q
                            cons[Locb]+=coe_phik[t][y]*G0[t][j][p,q]
                        else
                            cons[Locb]+=coe_phik[t][y]*2*G0[t][j][p,q]
                        end
                    end
                end
            end
        end


        for i=1:m[t]
            G[t][i]=Array{Any}(undef, lblock_sigma[t][i])
            for j=1:lblock_sigma[t][i]
                G[t][i][j]=@variable(model, [1:lt_block_sigma[t][i][j], 1:lt_block_sigma[t][i][j]],PSD)
                for p=1:lt_block_sigma[t][i][j]
                    for q=p:lt_block_sigma[t][i][j]
                        for z=1:lg[J[t][i]]
                            for y=1:lphik[t]
                                nota=supp_phik[t][:,y]+basis_sigma[t][i][:,block_sigma[t][i][j][p]]+basis_sigma[t][i][:,block_sigma[t][i][j][q]]+supp_g[J[t][i]][:,z]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if p==q
                                  cons[Locb]+=coe_phik[t][y]*coe_g[J[t][i]][z]*G[t][i][j][p,q]
                                else
                                  cons[Locb]+=coe_phik[t][y]*2*coe_g[J[t][i]][z]*G[t][i][j][p,q]
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
                H[t][i][j]=@variable(model, [1:lt_block_psi[t][i][j], 1:lt_block_psi[t][i][j]],Symmetric)
                for p=1:lt_block_psi[t][i][j]
                    for q=p:lt_block_psi[t][i][j]
                        for z=1:lh[W[t][i]]
                            for y=1:lphik[t]
                                nota=supp_phik[t][:,y]+basis_psi[t][i][:,block_psi[t][i][j][p]]+basis_psi[t][i][:,block_psi[t][i][j][q]]+supp_h[W[t][i]][:,z]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if p==q
                                  cons[Locb]+=coe_phik[t][y]*coe_h[W[t][i]][z]*H[t][i][j][p,q]
                                else
                                  cons[Locb]+=coe_phik[t][y]*2*coe_h[W[t][i]][z]*H[t][i][j][p,q]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    bc=[AffExpr(0) for i=1:llargest_supp]
    for i=1:lF
        Locb=bfind(largest_supp,llargest_supp,supp_F[:,i],n0)
        bc[Locb]=coe_F[i]
    end
    theta_bc=[AffExpr(0) for i=1:llargest_supp]
    for i=1:lphik0
        Locb=bfind(largest_supp,llargest_supp,supp_phik0[:,i],n0)
        theta_bc[Locb]=coe_phik0[i]
    end
    

    @variable(model, lambda)
    @constraint(model, bc-lambda.*theta_bc-cons.==0)
    @objective(model, Max, lambda)
    optimize!(model)

    opt_val = value(lambda)
    println("Termination status = ", termination_status(model))
    println("Optimal value = ",opt_val)


    #extraction of optimizers
    have_sol=1
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
        if sol[t]==[]
            have_sol=0
        end
    end
    
    if have_sol==0
        sol=[]
    end

    return opt_val, sol

end



function mix_sparsity_adding_spherical_constraints(x,g,h,I,J,W,eps,k,r)
    # small parameter
    p=length(I)
    n=Array{Int64}(undef, p)
    a0=Array{Any}(undef, p)
    a=Array{Any}(undef, p)
    omega0=Array{Any}(undef, p)
    omega=Array{Any}(undef, p)
    
    # Define centers and square of radius
    for t=1:p
        n[t]=length(I[t])
        a0[t]=zeros(Float32,n[t],1); a[t] = Matrix{Float32}(LinearAlgebra.I, n[t], n[t])
        omega0[t]=0; omega[t] = zeros(Float32,n[t], 1)
    end
    
    n_max=maximum(n)
    Lambda=Array{Any}(undef, n_max)
    for t=1:n_max
        Lambda[t]=0
        for j=1:p
            if n[j]>=t
                if Lambda[t]==0
                    Lambda[t]=[j]
                else
                    Lambda[t]=[Lambda[t];j]
                end
            end
        end
    end

    println("Determine L0:")
    println("------------------------------------")

    # Polynomial to optimize
    f=Array{Any}(undef,p,1)
    f[1]=sum(x[I[1]].^2)
    for t=2:p
        R=I[1]
        for j=2:t-1
            R=[R;I[j]]
        end
        R=setdiff(I[t], R)
        if length(R)>0
            f[t]=sum(x[R].^2)
        else
            f[t]=x[I[t][1]]-x[I[t][1]]
        end
    end
    
    f0 = sum(x.^2)

    L0,sol=PuVa_mix_hierarchy(x,f,f0,g,I,J,W,h,eps,k,r)
    
    println("------------------------------------")
    println("L0 = ", L0)
    println("====================================")
    
    if sol==[]
        for j=1:p
            #inequalities polynomial
            for i=1:p
                if j==1
                    if g==[]
                        g=[L0-sum(x[I[i]].^2)]
                    else
                        g=[g;L0-sum(x[I[i]].^2)]                       
                    end
                    if J[i]==[]
                        J[i]=[length(g)]
                    else
                        J[i]=[J[i];length(g)]
                    end  
                else
                    if i<j
                        g=g[setdiff(1:length(g), [length(g)-p+i])]
                        J[i]=J[i][1:end-p+i]
                        
                        if h==[]
                            h=[omega0[i]-sum((x[I[i]]-a0[i]).^2)]                        
                        else
                            h=[h;omega0[i]-sum((x[I[i]]-a0[i]).^2)]                       
                        end
                        if W[i]==[]
                            W[i]=[length(h)]
                        else
                            W[i]=[W[i];length(h)]
                        end   
                    else
                        J[i][j:end]=J[i][j:end].-1
                    end
                    
                end
            end
            f=Array{Any}(undef,p,1)
            for t=1:p
                if t!=j
                    f[t]=x[I[t][1]]-x[I[t][1]]
                else
                    f[t]=sum((x[I[j]]-a0[j]).^2)
                end
            end
            
            f0=f[j]
            println("Determine omega",0,"[",j,"]:")
            println("------------------------------------")
            omega0[j],sol=Pu_mix_hierarchy(x,f,f0,g,I,J,W,h,k,r)
            println("------------------------------------")
            println("omega",0,"[",j,"] = ", omega0[j])
            println("====================================")
            if sol!=[]
                break
            end              
        end



        #inequalities polynomial

        if sol==[]
            for t=1:n_max               
                for j in Lambda[t]
                    #inequalities polynomial
                    for i in Lambda[t]
                        if i<j
                            if h==[]
                                h=[omega[t][i]-sum((x[I[i]]-a[t][:,i]).^2)]                       
                            else
                                h=[h;omega[t][i]-sum((x[I[i]]-a[t][:,i]).^2)]                          
                            end
                            if W[i]==[]
                                W[i]=[length(h)]
                            else
                                W[i]=[W[i];length(h)]
                            end
                        end
                    end
                    f=Array{Any}(undef,p,1)
                    for t=1:p
                        if t!=j
                            f[t]=x[I[t][1]]-x[I[t][1]]
                        else
                            f[t]=sum((x[I[j]]-a0[j]).^2)
                        end
                    end
                    f0=f[j]
                    println("Determine omega",t,"[",j,"]:")
                    println("------------------------------------")

                    omega[t][j],sol=Pu_mix_hierarchy(x,f,f0,g,I,J,W,h,k,r)
                    println("------------------------------------")
                    println("omega",t,"[",j,"] = ", omega[t][j])
                    println("====================================")
                    if sol!=[]
                        break
                    end

                end
                if sol!=[]
                    break
                end            
            end
        end
    end
    return sol
end