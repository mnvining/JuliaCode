using SymPy
function hsol(z,x,Al)
    # z is grid of pts
    x=Sym("x")
    f(x)=exp(1/sqrt(Al)*(1-x))
    g(x)=exp(1/sqrt(Al)*(x-1))
        return h_sol_1(x)=exp(1/sqrt(Al)*(1-x))/norm(evalf(f(x),z)), h_sol_2(x)=exp(1/sqrt(Al)*(x-1))/norm(evalf(g(x),z))
     
end
