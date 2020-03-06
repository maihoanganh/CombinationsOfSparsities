function Pu_sparse_hierarchy(x,f0,g,h,I,J,W,k)

    n0=length(x) # Number of variables
    m0=length(g) # Number of constraints
    l0=length(h) # Number of constraints

    lf0,supp_f0,coe_f0=sparse_info(f0,x,n0)

    supp_g=Vector{SparseMatrixCSC}(undef, m0)
    coe_g=Vector{Vector{Float64}}(undef, m0)
    lg=Vector{UInt16}(undef, m0)
    for j=1:m0
        lg[j],supp_g[j],coe_g[j]=sparse_info(g[j],x,n0)
    end

    supp_h=Vector{SparseMatrixCSC}(undef, l0)
    coe_h=Vector{Vector{Float64}}(undef, l0)
    lh=Vector{UInt16}(undef, l0)
    for j=1:l0
        lh[j],supp_h[j],coe_h[j]=sparse_info(h[j],x,n0)
    end

    p=length(I)

    n=Vector{Int64}(undef, p)
    m=Vector{UInt16}(undef, p)
    l=Vector{UInt16}(undef, p)

    u=Vector{Vector{UInt16}}(undef, p)
    w=Vector{Vector{UInt16}}(undef, p)

    lu0=Vector{UInt32}(undef, p)
    lu=Vector{Vector{UInt32}}(undef, p)
    lw=Vector{Vector{UInt32}}(undef, p)

    basis_sigma0=Vector{SparseMatrixCSC}(undef, p)
    basis_sigma=Vector{Vector{SparseMatrixCSC}}(undef, p)
    basis_psi=Vector{Vector{SparseMatrixCSC}}(undef, p)

    largest_supp=sparse(UInt16[])

    #Degrees of constraints
    for t=1:p
        n[t]=length(I[t])

        lu0[t]=binomial(n[t]+k,n[t])
        basis_sigma0[t]=spzeros(n0,lu0[t])
        basis_sigma0[t][I[t],:]=get_basis(n[t],k)
        

        for i=1:lu0[t]
            for r=i:lu0[t]
                nota=basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                if largest_supp==[]
                    @inbounds largest_supp=nota
                else
                    @inbounds largest_supp=[largest_supp nota]
                end
            end
        end

        m[t]=length(J[t])
        u[t]=Vector{UInt16}(undef, m[t])
        lu[t]=Vector{UInt32}(undef, m[t])

        basis_sigma[t]=Vector{SparseMatrixCSC}(undef, m[t])
        for j=1:m[t]
            u[t][j]=ceil(UInt16,.5*maxdegree(g[J[t][j]]))
            lu[t][j]=binomial(n[t]+k-u[t][j],n[t])
            @inbounds basis_sigma[t][j]=basis_sigma0[t][:,1:lu[t][j]]

            for i=1:lu[t][j]
                for r=i:lu[t][j]
                    for z=1:lg[J[t][j]]
                        @inbounds nota=supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                        @inbounds largest_supp=[largest_supp nota]
                    end
                end
            end
        end

        l[t]=length(W[t])
        w[t]=Vector{UInt16}(undef, l[t])
        lw[t]=Vector{UInt32}(undef, l[t])

        basis_psi[t]=Vector{SparseMatrixCSC}(undef, l[t])
        for j=1:l[t]
            w[t][j]=ceil(UInt16,.5*maxdegree(h[W[t][j]]))
            lw[t][j]=binomial(n[t]+k-w[t][j],n[t])
            @inbounds basis_psi[t][j]=basis_sigma0[t][:,1:lw[t][j]]

            for i=1:lw[t][j]
                for r=i:lw[t][j]
                    for z=1:lh[W[t][j]]
                        @inbounds nota=supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                        @inbounds largest_supp=[largest_supp nota]
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

    G0=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, p)
    G=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, p)
    H=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, p)

    for t=1:p
        if lu0[t]==1
            @inbounds G0[t]=@variable(model, lower_bound=0)
            @inbounds nota=UInt8(2)*basis_sigma0[t]
            Locb=bfind(largest_supp,llargest_supp,nota,n0)
            @inbounds add_to_expression!(cons[Locb],G0[t])
        else
            G0[t]=@variable(model, [1:lu0[t], 1:lu0[t]],PSD)
            for i=1:lu0[t]
                for r=i:lu0[t]
                    @inbounds nota=basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                    Locb=bfind(largest_supp,llargest_supp,nota,n0)
                    if i==r
                      @inbounds add_to_expression!(cons[Locb],G0[t][i,r])
                    else
                      @inbounds add_to_expression!(cons[Locb],2*G0[t][i,r])
                    end
                end
            end
        end

        G[t]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, m[t])
        for j=1:m[t]
            if lu[t][j]==1
                @inbounds G[t][j]=@variable(model, lower_bound=0)
                for z=1:lg[J[t][j]]
                    @inbounds nota=supp_g[J[t][j]][:,z]+UInt8(2)*basis_sigma[t][j]
                    Locb=bfind(largest_supp,llargest_supp,nota,n0)
                    @inbounds add_to_expression!(cons[Locb],G[t][j]*coe_g[J[t][j]][z])
                end
            else
                G[t][j]=@variable(model, [1:lu[t][j], 1:lu[t][j]],PSD)
                for i=1:lu[t][j]
                    for r=i:lu[t][j]
                        for z=1:lg[J[t][j]]
                            @inbounds nota=supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                            Locb=bfind(largest_supp,llargest_supp,nota,n0)
                            if i==r
                              @inbounds add_to_expression!(cons[Locb],coe_g[J[t][j]][z]*G[t][j][i,r])
                            else
                              @inbounds add_to_expression!(cons[Locb],2*coe_g[J[t][j]][z]*G[t][j][i,r])
                            end
                        end
                    end
                end
            end
        end

        H[t]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, l[t])
        for j=1:l[t]
           if lw[t][j]==1
                @inbounds H[t][j]=@variable(model)
                for z=1:lh[W[t][j]]
                    @inbounds nota=supp_h[W[t][j]][:,z]+UInt8(2)*basis_psi[t][j]
                    Locb=bfind(largest_supp,llargest_supp,nota,n0)
                    @inbounds add_to_expression!(cons[Locb],H[t][j]*coe_h[W[t][j]][z])
                end
            else
                H[t][j]=@variable(model, [1:lw[t][j], 1:lw[t][j]],Symmetric)
                for i=1:lw[t][j]
                    for r=i:lw[t][j]
                        for z=1:lh[W[t][j]]
                            @inbounds nota=supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                            Locb=bfind(largest_supp,llargest_supp,nota,n0)
                            if i==r
                              @inbounds add_to_expression!(cons[Locb],coe_h[W[t][j]][z]*H[t][j][i,r])
                            else
                              @inbounds add_to_expression!(cons[Locb],2*coe_h[W[t][j]][z]*H[t][j][i,r])
                            end
                        end
                    end
                end
            end
        end
    end

    bc=[AffExpr(0) for i=1:llargest_supp]
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
    flag=1
    sol=Vector{Vector{Float64}}(undef,p)
    for t=1:p
        println("====================================")
        println("Atom(I[",t,"]):")
        println("====================================")
        sol[t]=extract_optimizers(value.(G0[t]),lu0[t],basis_sigma0[t],n0,n[t],m[t],l[t],opt_val,0,g[J[t]],h[W[t]],x,I[t])
        if sol[t]==[]
            flag=0
        end
    end
    
    if flag==0
        sol=Vector{Float64}[]
    end

    return opt_val,sol

end





function PuVa_sparse_hierarchy(x,f0,g,h,I,J,W,eps,k)

    n0=length(x) # Number of variables
    m0=length(g) # Number of constraints
    l0=length(h) # Number of constraints
    
    d=1+floor(Int16,.5*maxdegree(f0))
    p=length(I)
    
    #Determine hat I_j and s_j
    hat_I=Vector{Vector{UInt8}}(undef, p)
    s=Vector{UInt8}(undef, p)
    theta=Vector{Polynomial{true,Int16}}(undef, p)
    hat_theta=Vector{Polynomial{true,Int16}}(undef, p)
    for j=1:p
        if j>=2
            union_I=UInt8(0)
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
    D=Vector{Polynomial{true,Int16}}(undef, p)
    Theta=Vector{Polynomial{true,Int16}}(undef, p)
    omega=Vector{UInt16}(undef, p)
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
        omega[j]=ceil(UInt16,.5*maxdegree(Theta[j]))
    end
    
    

    psi=0
    phi0=1
    phi=Vector{Polynomial{true,Int16}}(undef, p)
    
    supp_phik=Vector{SparseMatrixCSC}(undef, p)
    coe_phik=Vector{Vector{Float64}}(undef, p)
    lphik=Vector{UInt16}(undef, p)
    for j=1:p
        psi+=theta[j]^d
        phi0 *=Theta[j]

        phi[j]=1
        for r=1:p
            if j!=r
                phi[j]*=Theta[r]
            end
        end
        lphik[j],supp_phik[j],coe_phik[j]=sparse_info(phi[j]^k,x,n0)
    end
    ############
    lphik0,supp_phik0,coe_phik0=sparse_info(phi0^k,x,n0)
   

    F=phi0^k*(f0+eps*psi)
    
    lF,supp_F,coe_F=sparse_info(F,x,n0)

    supp_g=Vector{SparseMatrixCSC}(undef, m0)
    coe_g=Vector{Vector{Float64}}(undef, m0)
    lg=Vector{UInt16}(undef, m0)
    for j=1:m0
        lg[j],supp_g[j],coe_g[j]=sparse_info(g[j],x,n0)
    end
    ###############
    
    
    supp_h=Vector{SparseMatrixCSC}(undef, l0)
    coe_h=Vector{Vector{Float64}}(undef, l0)
    lh=Vector{UInt16}(undef, l0)
    for j=1:l0
        lh[j],supp_h[j],coe_h[j]=sparse_info(h[j],x,n0)
    end

    
    n=Vector{Int64}(undef, p)
    m=Vector{UInt16}(undef, p)
    l=Vector{UInt16}(undef, p)

    u=Vector{Vector{UInt16}}(undef, p)
    w=Vector{Vector{UInt16}}(undef, p)

    lu0=Vector{UInt32}(undef, p)
    lu=Vector{Vector{UInt32}}(undef, p)
    lw=Vector{Vector{UInt32}}(undef, p)

    basis_sigma0=Vector{SparseMatrixCSC}(undef, p)
    basis_sigma=Vector{Vector{SparseMatrixCSC}}(undef, p)
    basis_psi=Vector{Vector{SparseMatrixCSC}}(undef, p)

    largest_supp=sparse(UInt16[])
    
    #Degrees of constraints
    for t=1:p
        n[t]=length(I[t])

        lu0[t]=binomial(n[t]+d+k*omega[t],n[t])
        basis_sigma0[t]=spzeros(n0,lu0[t])
        basis_sigma0[t][I[t],:]=get_basis(n[t],d+k*omega[t])

        for i=1:lu0[t]
            for r=i:lu0[t]
                for z=1:lphik[t]
                    nota=supp_phik[t][:,z]+basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                    if largest_supp==[]
                        largest_supp=nota
                    else
                        largest_supp=[largest_supp nota]
                    end
                end
            end
        end

        m[t]=length(J[t])
        u[t]=Vector{UInt16}(undef, m[t])
        lu[t]=Vector{UInt32}(undef, m[t])

        basis_sigma[t]=Vector{SparseMatrixCSC}(undef, m[t])
        for j=1:m[t]
            u[t][j]=ceil(UInt16,.5*maxdegree(g[J[t][j]]))
            lu[t][j]=binomial(n[t]+d+k*omega[t]-u[t][j],n[t])
            basis_sigma[t][j]=basis_sigma0[t][:,1:lu[t][j]]

            for i=1:lu[t][j]
                for r=i:lu[t][j]
                    for z=1:lg[J[t][j]]
                        for q=1:lphik[t]
                            nota=supp_phik[t][:,q]+supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                            largest_supp=[largest_supp nota]
                        end
                    end
                end
            end
        end

        l[t]=length(W[t])
        w[t]=Vector{UInt16}(undef, l[t])
        lw[t]=Vector{UInt32}(undef, l[t])

        basis_psi[t]=Vector{SparseMatrixCSC}(undef, l[t])
        for j=1:l[t]
            w[t][j]=ceil(UInt16,.5*maxdegree(h[W[t][j]]))
            lw[t][j]=binomial(n[t]+d+k*omega[t]-w[t][j],n[t])
            basis_psi[t][j]=basis_sigma0[t][:,1:lw[t][j]]

            for i=1:lw[t][j]
                for r=i:lw[t][j]
                    for z=1:lh[W[t][j]]
                        for q=1:lphik[t]
                            nota=supp_phik[t][:,q]+supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                            largest_supp=[largest_supp nota]
                        end
                    end
                end
            end
        end

    end
    ##############
    

    largest_supp=unique(largest_supp,dims=2)
    largest_supp=sortslices(largest_supp,dims=2)
    llargest_supp=size(largest_supp,2)

    model=Model(with_optimizer(Mosek.Optimizer, QUIET=true))
    cons=[AffExpr(0) for i=1:llargest_supp]

    G0=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, p)
    G=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, p)
    H=Vector{Vector{Union{VariableRef,Symmetric{VariableRef}}}}(undef, p)

    for t=1:p
        if lu0[t]==1
            @inbounds G0[t]=@variable(model, lower_bound=0)
            for q=1:lphik[t]
                @inbounds nota=supp_phik[t][:,q]+UInt8(2)*basis_sigma0[t]
                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*G0[t])
            end
        else
            G0[t]=@variable(model, [1:lu0[t], 1:lu0[t]],PSD)
            for i=1:lu0[t]
                for r=i:lu0[t]
                    for q=1:lphik[t]
                        @inbounds nota=supp_phik[t][:,q]+basis_sigma0[t][:,i]+basis_sigma0[t][:,r]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        if i==r
                          @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*G0[t][i,r])
                        else
                          @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*2*G0[t][i,r])
                        end
                    end
                end
            end
        end

        G[t]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, m[t])
        for j=1:m[t]
            if lu[t][j]==1
                @inbounds G[t][j]=@variable(model, lower_bound=0)
                for z=1:lg[J[t][j]]
                    for q=1:lphik[t]
                        @inbounds nota=supp_phik[t][:,q]+supp_g[J[t][j]][:,z]+UInt8(2)*basis_sigma[t][j]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*G[t][j])
                    end
                end
            else
                G[t][j]=@variable(model, [1:lu[t][j], 1:lu[t][j]],PSD)
                for i=1:lu[t][j]
                    for r=i:lu[t][j]
                        for z=1:lg[J[t][j]]
                            for q=1:lphik[t]
                                nota=supp_phik[t][:,q]+supp_g[J[t][j]][:,z]+basis_sigma[t][j][:,i]+basis_sigma[t][j][:,r]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if i==r
                                  @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*coe_g[J[t][j]][z]*G[t][j][i,r])
                                else
                                  @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*2*coe_g[J[t][j]][z]*G[t][j][i,r])
                                end
                            end
                        end
                    end
                end
            end
        end

        H[t]=Vector{Union{VariableRef,Symmetric{VariableRef}}}(undef, l[t])
        for j=1:l[t]
            if lw[t][j]==1
                @inbounds G[t][j]=@variable(model)
                for z=1:lh[W[t][j]]
                    for q=1:lphik[t]
                        @inbounds nota=supp_phik[t][:,q]+supp_h[W[t][j]][:,z]+UInt8(2)*basis_psi[t][j]
                        Locb=bfind(largest_supp,llargest_supp,nota,n0)
                        @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*H[t][j])
                    end
                end
            else
                H[t][j]=@variable(model, [1:lw[t][j], 1:lw[t][j]],Symmetric)
                for i=1:lw[t][j]
                    for r=i:lw[t][j]
                        for z=1:lh[W[t][j]]
                            for q=1:lphik[t]
                                @inbounds nota=supp_phik[t][:,q]+supp_h[W[t][j]][:,z]+basis_psi[t][j][:,i]+basis_psi[t][j][:,r]
                                Locb=bfind(largest_supp,llargest_supp,nota,n0)
                                if i==r
                                  @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*coe_h[W[t][j]][z]*H[t][j][i,r])
                                else
                                  @inbounds add_to_expression!(cons[Locb],coe_phik[t][q]*2*coe_h[W[t][j]][z]*H[t][j][i,r])
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    ##################
    
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
    flag=1
    sol=Vector{Vector{Float64}}(undef,p)
    for t=1:p
        println("====================================")
        println("Atom(I[",t,"]):")
        println("====================================")
        sol[t]=extract_optimizers(value.(G0[t]),lu0[t],basis_sigma0[t],n0,n[t],m[t],l[t],opt_val,0,g[J[t]],h[W[t]],x,I[t])
        if sol[t]==[]
            flag=0
        end
    end
    
    if flag==0
        sol=Vector{Float64}[]
    end

    return opt_val,sol
end



function sparse_ASC(x,g,h,I,J,W,eps,k)
    # small parameter
    p=length(I)
    n=Vector{UInt32}(undef, p)
    a=Vector{Matrix{Float32}}(undef, p)
    omega0=Vector{Float32}(undef, p)
    omega=Vector{Vector{Float32}}(undef, p)
    
    # Define centers and square of radius
    for t=1:p
        n[t]=length(I[t])
        a[t] = Matrix{Float32}(LinearAlgebra.I, n[t], n[t])
        omega0[t]=0; omega[t] = zeros(Float32,n[t])
    end
    
    n_max=maximum(n)
    Lambda=Vector{Vector{UInt32}}(undef, n_max)
    for t=1:n_max
        Lambda[t]=[0]
        for j=1:p
            if n[j]>=t
                if Lambda[t]==[0]
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
    f0 = sum(x.^2)

    L0,sol=PuVa_sparse_hierarchy(x,f0,g,h,I,J,W,eps,k)
    
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
                            h=[omega0[i]-sum(x[I[i]].^2)]                        
                        else
                            h=[h;omega0[i]-sum(x[I[i]].^2)]                       
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
            f0=sum(x[I[j]].^2)
            println("Determine omega",0,"[",j,"]:")
            println("------------------------------------")
            omega0[j],sol=Pu_sparse_hierarchy(x,f0,g,h,I,J,W,k)
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
                    f0=sum((x[I[j]]-a[t][:,j]).^2)
                    println("Determine omega",t,"[",j,"]:")
                    println("------------------------------------")

                    omega[t][j],sol=Pu_sparse_hierarchy(x,f0,g,h,I,J,W,k)
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


function extract_infinite_set_optimizers(x,lambda,f,g,I,J,W,h,eps,k)
    n=length(x)
    p=length(I)
    @polyvar y[1:p-1,1]
    X=[x;y]
    
    I[1]=[I[1];n+1]
    I[p]=[I[1];n+p-1]
    
    g=[g;lambda-f[p]-X[n+p-1]]
    J[p]=[J[p];length(g)]
    
    h=[h;f[1]-X[n+1]]
    W[1]=[W[1];length(h)]
    
    for j=2:p-1
        I[j]=[I[j];n+j-1;n+j]
        
        h=[h;X[n+j-1]+f[j]-X[n+j]]
        W[j]=[W[j];length(h)]
    end
    
    sol=sparse_ASC(X,g,h,I,J,W,eps,k)
    return sol
end
    