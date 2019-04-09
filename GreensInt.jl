    using SymPy
using GenericSVD
include("coeffcalc.jl")
function GreensInt(d,Al,p)
    x=symbols("x")
    if d==0
        return (-exp(-(x+1)/sqrt(Al))-exp(x/sqrt(Al))+exp(-1/sqrt(Al))+1)/(exp(-1/sqrt(Al))+1)
    elseif d==1
        return((exp(-1/sqrt(Al)))^2*x+exp(-1/sqrt(Al))*exp(x/sqrt(Al))-exp(-(x+1)/sqrt(Al))-x)/(-1+(exp(-1/sqrt(Al)))^2)
    elseif d==2
        return ((2*Al+1)*exp(-(x+1)/sqrt(Al))-2*Al*exp(-(x+2)/sqrt(Al))+(-1-2*Al)*exp((-1+x)/sqrt(Al))+(x^2+2*Al)*exp(-2/sqrt(Al))-x^2+2*Al*exp(x/sqrt(Al))-2*Al)/(-1+exp(-2/sqrt(Al)))
    elseif d==3
        return ((-6*Al-1)*exp(-(x+1)/sqrt(Al))+(6*Al+1)*exp((x-1)/sqrt(Al))+x*(-1+exp(-2/sqrt(Al)))*(x^2+6*Al))/(-1+exp(-2/sqrt(Al)))
    elseif d==4
        return ((24*Al^2+12*Al+1)*exp(-(x+1)/sqrt(Al))-24*Al^2*exp(-(2+x)/sqrt(Al))+(-24*Al^2-12*Al-1)*exp((x-1)/sqrt(Al))+(x^4+12*Al*x^2+24*Al^2)*exp(-2/sqrt(Al))-x^4-12*x^2*Al+24*Al^2*exp(x/sqrt(Al))-24*Al^2)/(-1+exp(-2/sqrt(Al)))
    elseif d==5
        return ((-120*Al^2-20*Al-1)*exp(-(x+1)/sqrt(Al))+(120*Al^2+20*Al+1)*exp((x-1)/sqrt(Al))+x*(x^4+20*Al*x^2+120*Al^2)*(-1+exp(-2/sqrt(Al))))/(-1+exp(-2/sqrt(Al)))
    elseif d==6
        return -(720*exp(-1/sqrt(Al))*Al^3-720*Al^3-360*Al^2-30*Al-1)*exp((-x-1)/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)-(720*exp(-1/sqrt(Al))*Al^3-720*Al^3+360*exp(-1/sqrt(Al))*Al^2+30*exp(-1/sqrt(Al))*Al+exp(-1/sqrt(Al)))*exp(x/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)+x^6+30*Al*x^4+360*Al^2*x^2+720*Al^3
    elseif d==7
        return ((-5040*Al^3-840*Al^2-42*Al-1)*exp(-(x+1)/sqrt(Al))+(5040*Al^3+840*Al^2+42*Al+1)*exp((x-1)/sqrt(Al))+x*(x^6+42*Al*x^4+840*Al^2*x^2+5040*Al^3)*(-1+exp(-2/sqrt(Al))))/(-1+exp(-2/sqrt(Al)))  
    elseif d==8
        return -(40320*exp(-1/sqrt(Al))*Al^4-40320*Al^4-20160*Al^3-1680*Al^2-56*Al-1)*exp((-x-1)/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)-(40320*exp(-1/sqrt(Al))*Al^4-40320*Al^4+20160*exp(-1/sqrt(Al))*Al^3+1680*exp(-1/sqrt(Al))*Al^2+56*exp(-1/sqrt(Al))*Al+exp(-1/sqrt(Al)))*exp(x/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)+x^8+56*Al*x^6+1680*Al^2*x^4+20160*Al^3*x^2+40320*Al^4
    elseif d==9        
        return (-362880*Al^4-60480*Al^3-3024*Al^2-72*Al-1)*exp(-(x+1)/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)+exp(-1/sqrt(Al))*(362880*Al^4+60480*Al^3+3024*Al^2+72*Al+1)*exp(x/sqrt(Al))/(-1+(exp(-1/sqrt(Al)))^2)+x^9+72*Al*x^7+3024*Al^2*x^5+60480*Al^3*x^3+362880*x*Al^4   
    end
end

function totalGreen(d,n,Al,p)
    x=symbols("x")
    expr=0*x
    zz=coeffcalc(d,n)
    for j=1:length(zz)
        m(x)=GreensInt(j-1,Al,p)
        expr=expr+zz[j]*m(x)
    end
    return expr
end

