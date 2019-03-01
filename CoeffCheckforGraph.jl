using PyPlot
using SymPy
using GenericSVD
using LinearAlgebra
include("All.jl")
include("dftmtx.jl")
include("coeffcalc.jl")
include("hsol.jl")
include("ConvertToBig.jl")
include("GreensInt.jl")
include("Spaces.jl")
include("RUNGraph.jl")
mm = 60
Al=logspace(-5,-1,mm)
StabMat=zeros(mm,1)
A1=zeros(mm,1)
A2=zeros(mm,1)
AccMat=zeros(mm,1)
Stab2Mat=zeros(mm,1)


for i=1:mm
    (a,b,c,d,f)=RUNGraph(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(2),1)
    StabMat[i]=a
    A1[i]=b
    A2[i]=c
    AccMat[i]=d
    Stab2Mat[i]=f
end
figure(1)
loglog(Al,abs.(A1),label="A1")
loglog(Al,abs.(A2),label="A2")
loglog(Al,StabMat,label="Stability")
loglog(Al,Stab2Mat,label="Greens")
loglog(Al,ones(size(Al)),label="f(a)=1")
       legend()
title("Function")

for i=1:mm
    (a,b,c,d,f)=RUNGraph(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(2),2)
    StabMat[i]=a
    A1[i]=b
    A2[i]=c
    AccMat[i]=d
    Stab2Mat[i]=f
end
figure(2)
loglog(Al,abs.(A1),label="A1")
loglog(Al,abs.(A2),label="A2")
loglog(Al,StabMat,label="Stability")
loglog(Al,Stab2Mat,label="Greens")
loglog(Al,AccMat,label="Accuracy")
loglog(Al,ones(size(Al)),label="f(a)=1")
              legend()
title("2nd Derivative")

figure(3)
for i=1:mm
    (a,b,c,d,f)=RUNGraph(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(2),3)
    StabMat[i]=a
    A1[i]=b
    A2[i]=c
    AccMat[i]=d
    Stab2Mat[i]=f
end
loglog(Al,abs.(A1),label="A1")
loglog(Al,abs.(A2),label="A2")
loglog(Al,StabMat,label="Stability")
loglog(Al,Stab2Mat,label="Greens")
loglog(Al,ones(size(Al)),label="f(a)=1")
legend()
title("Modified Sobolev")



