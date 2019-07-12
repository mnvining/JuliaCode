

include("MP.jl")
include("IDOCalc.jl")
include("NCalc.jl")
using PyPlot
Al=logspace(-8,0,27)
w=logspace(0,9,36)
yc=zeros(31,36)
ys=zeros(30,36)
Gc=zeros(31,36)
Gs=zeros(30,36)
cosacc=zeros(27,36)
sinacc=zeros(27,36)
cosstab=zeros(27,36)
sinstab=zeros(27,36)

for i=1:27
    (IDOC,IDOS)=IDOCalc(Al[i])
    for j=1:36
        (yc[:,j],ys[:,j],Gc[:,j],Gs[:,j])=MP(0,w[j],FCMat,FCMatS,C_C,C_S,IDOC,IDOS)
        (cosacc[i,j],sinacc[i,j],cosstab[i,j],sinstab[i,j])=NCalc(0,w[j],FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm)
    end
    println(i)
end
