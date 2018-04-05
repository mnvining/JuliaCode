    using SymPy

include("coeffcalc.jl")
include("All.jl")
function Greens(x,a,Al)


    A1=-(1/2)*exp(-4/sqrt(Al))*sqrt(Al)*(-exp((a-1)/sqrt(Al))+exp(-(a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1)
    A2=(1/2)*sqrt(Al)*(-exp((a-1)/sqrt(Al))+exp(-(a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1)
    B1=-(1/2)*sqrt(Al)*(exp(-(a-1)/sqrt(Al))*exp(-4/sqrt(Al))-exp((a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1)
    B2=(1/2)*sqrt(Al)*(exp(-(a-1)/sqrt(Al))*exp(-4/sqrt(Al))-exp((a-1)/sqrt(Al)))/(exp(-4/sqrt(Al))-1)

        return piecewise((A1*exp((1-x)/sqrt(Al))+A2*exp((x-1)/sqrt(Al)),(Lt(x,a))),(B1*exp((1-x)/sqrt(Al))+B2*exp((x-1)/sqrt(Al)),(Ge(x,a))))

end

function CalcGreens(d,n,Al)
    x,a=Sym("x","a")
    G(x,a)=Greens(x,a,Al)
    f(x)=fmaker(d,n) 
    return integrate(f(a)*G(x,a),(a,-1,1))#+integrate(f(a)*G(x,a),(a,x,1))
end









