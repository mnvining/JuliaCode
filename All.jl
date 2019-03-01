     using SymPy
using LinearAlgebra


function fmaker(d,n)
    x=Sym("x")
    h=2/(n-1);
    z=collect(-1:h:1-h);
    A=fliplr(myvander(z))
    Q,R=qr(A)
    if d==0
        expr=1/R[1,1]*x^0
        return expr
    else
        expr=x^d
        for i=0:d-1
            qi=fmaker(i,n)
            expr=expr-qi*R[i+1,d+1]
        end
        return expr/(R[d+1,d+1])
    end
end

    

      
