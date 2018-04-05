using SymPy
function hsol(z,Al)
    # z is grid of pts
        x = Sym("x")
    function h_sol_1(x)
  
      return exp(1/sqrt(Al)*(1-x))/norm(evalf(exp(1/sqrt(Al)*(1-x)),z))
    end

    function h_sol_2(x)

      return exp(1/sqrt(Al)*(x-1))/norm(evalf(exp(1/sqrt(Al)*(x-1)),z))
    end

    return h_sol_1, h_sol_2
end
