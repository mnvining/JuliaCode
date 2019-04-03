using SymPy
using LinearAlgebra
using PyPlot
using GenericSVD

function StabBisect(d,Al)
    include("CosCont.jl")
    Tol=Al^4;
    V=10
    a=BigFloat(-10)
    b=BigFloat(-7)
    for i=1:100
        Q=(b+a)/2;
        (a1,s1)=CosCont(d,Al,BigFloat(10)^a);
        (a2,s2)=CosCont(d,Al,BigFloat(10)^b);
        (a3,s3)=CosCont(d,Al,BigFloat(10)^Q);

        if s3 > 1-Tol
            if s1 < 1-Tol
                b=Q
            else
                a=Q
            end
            
        else
            if s2 < 1 - Tol
                b=Q
            else
                a=Q
            end   
        end
        println(Q)
        V=1-s2
        if V<Tol
            return Q
            break
            break
        end
        
        
    end
   

end

    
        
