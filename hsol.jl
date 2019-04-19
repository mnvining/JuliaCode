using SymPy
function hsol(z,Al,p)
    # z is grid of pts
    x=Sym("x")

    h_sol_1(x)=exp(1/sqrt(Al)*(-x-1))
    h_sol_2(x)=exp(1/sqrt(Al)*(x-p))

    return h_sol_1,h_sol_2
end

function evalasarray(expr,z)
    # expr is expression in symbolic using SymPy, x
    # z is array
    # N only works on single values, this returns array of evaluations
    x=Sym("x")
    q=zeros(BigFloat,size(z));
    for i=1:length(z)
        q[i]=subs(expr,x,z[i])
    end
    return q
end
