function info(f,x,n)
    mon=monomials(f)
    coe=coefficients(f)
    lmon=length(mon)
    supp=zeros(UInt8,n,lmon)
    for i=1:lmon
        for j in 1:n
            @inbounds supp[j,i]=MultivariatePolynomials.degree(mon[i],variable(x[j]))
        end
    end
    return lmon, supp, coe
end
#######################################

function sparse_info(f,x,n)
    mon=monomials(f)
    coe=coefficients(f)
    lmon=length(mon)
    supp=spzeros(UInt8,n,lmon)
    for i=1:lmon
        for j in 1:n
            @inbounds supp[j,i]=MultivariatePolynomials.degree(mon[i],variable(x[j]))
        end
    end
    return lmon, supp, coe
end

function get_basis(n,d)
    lb=binomial(n+d,d)
    basis=zeros(UInt8,n,lb)
    i=0
    t=1
    while i<d+1
        if basis[n,t]==i
           if i<d
              t=t+1
              basis[1,t]=i+1
              i=i+1
           else i=i+1
           end
        else j=1
             while basis[j,t]==0
                   j=j+1
             end
             if j==1
                t=t+1
                basis[:,t]=basis[:,t-1]
                basis[1,t]=basis[1,t]-1
                basis[2,t]=basis[2,t]+1
             else t=t+1
                  basis[:,t]=basis[:,t-1]
                  basis[1,t]=basis[j,t]-1
                  basis[j,t]=0
                  basis[j+1,t]=basis[j+1,t]+1
             end
        end
    end
    return basis
end


function bfind(A,l,a,n)
    if l==0
        return 0
    end
    low=1
    high=l
    while low<=high
        mid=Int(ceil(1/2*(low+high)))
        order=comp(A[:,mid],a,n)
        if order==0
           return mid
        elseif order<0
           low=mid+1
        else
           high=mid-1
        end
    end
    return 0
end

function comp(a,b,n)
    i=1
    while i<=n
          if a[i]<b[i]
             return -1
          elseif a[i]>b[i]
             return 1
          else
             i+=1
          end
    end
    if i==n+1
       return 0
    end
end

