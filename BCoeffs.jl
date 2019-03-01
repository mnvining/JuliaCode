using SymPy
using LinearAlgebra
using GenericSVD
using PyPlot
function BCoeffCalc(Fine_Grid,Coarse_Grid,Al,d,opt)
# opt: 1 for func
# opt: 2 for 2nd der
# opt: 3 for Sobolev-ish
include("GreensInt.jl")
include("Spaces.jl")
include("hsol.jl")
include("All.jl")
include("ConvertToBig.jl")
x=Sym("x")
n=length(Fine_Grid)

GG(x)=totalGreen(d,n,Al)
    (h1,h2)=hsol(Coarse_Grid,Al)
    h1Prime(x)=diff(h1(x),x,1)
    h2Prime(x)=diff(h2(x),x,1)
    d2h1(x)=diff(h1(x),x,2)
    d2h2(x)=diff(h2(x),x,2)
    d2GG(x)=diff(GG(x),x,2)

    
    EGG=ConvertToBig(evalasarray(GG(x),Coarse_Grid))
    Eh1=ConvertToBig(evalasarray(h1(x),Coarse_Grid))
    Eh2=ConvertToBig(evalasarray(h2(x),Coarse_Grid))

    EG2=ConvertToBig(evalasarray(d2GG(x),Coarse_Grid))
    EH1=ConvertToBig(evalasarray(d2h1(x),Coarse_Grid))
    EH2=ConvertToBig(evalasarray(d2h2(x),Coarse_Grid))

    
    EGGf=ConvertToBig(evalasarray(GG(x),Fine_Grid))
    Eh1f=ConvertToBig(evalasarray(h1(x),Fine_Grid))
    Eh2f=ConvertToBig(evalasarray(h2(x),Fine_Grid))

if opt==1
    Z=[Eh1'*Eh1 Eh2'*Eh1; Eh1'*Eh2 Eh2'*Eh2]
    zz=[-EGG'*Eh1;-EGG'*Eh2]
    cobb=\(Z,zz)

elseif opt==2
    Z=[EH1'*EH1 EH2'*EH1; EH1'*EH2 EH2'*EH2]
    zz=[-EG2'*EH1;-EG2'*EH2]
    cobb=\(Z,zz)

    elseif opt==3
    Z=[Eh1'*Eh1+Eh1'*EH1+EH1'*Eh1+EH1'*EH1 Eh2'*Eh1+EH2'*Eh1+Eh2'*EH1+EH2'*EH1; Eh1'*Eh2+Eh1'*EH2+EH1'*Eh2+EH1'*EH2 Eh2'*Eh2+EH2'*Eh2+Eh2'*EH2+EH2'*EH2]
    zz=[-(EGG'*Eh1+EG2'*Eh1+EGG'*EH1+EG2'*EH1);-(EGG'*Eh2+EG2'*Eh2+EGG'*EH2+EG2'*EH2)]
    cobb=\(Z,zz)
end
    return cobb,EGGf,Eh1f,Eh2f
end

