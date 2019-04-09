using SymPy
include("All.jl")

function coeffcalc(d,n)
    x=symbols("x")
    exx=fmaker(d,n)
    co=N(exx(x=>0))
    for i=1:d
        co=vcat(co,N(exx.coeff(x^i)))
    end
    return co
end

 
function fliplr(A)
    (m,n)=size(A)
    N=zeros(BigFloat,size(A))

    for j=1:n
        N[:,n-j+1]=A[:,j]
    end
    return N
end

function myvander(x)
    
    L=length(x)
    M=zeros(BigFloat,L,L)
    for i=1:L
        M[:,L-i+1]=x.^(i-1)
    end
    return M
end

    


      

        
