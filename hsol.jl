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
    q=zeros(size(z));
    for i=1:length(z)
        p=subs(expr,x,z[i])
        q[i]=N(p)
    end
    return q
end

        
