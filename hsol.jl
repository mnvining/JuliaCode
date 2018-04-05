using SymPy
function hsol(z,Al)
    # z is grid of pts
    x=Sym("x")
    Ar1=evalasarray(exp(1/sqrt(Al)*(1-x)),z)
    Ar2=evalasarray(exp(1/sqrt(Al)*(x-1)),z)
    val1=norm(Ar1)
    val2=norm(Ar2)
    h_sol_1(x)=exp(1/sqrt(Al)*(1-x))/val1
    h_sol_2(x)=exp(1/sqrt(Al)*(x-1))/val2

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

        
