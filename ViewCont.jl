

include("MP.jl")
include("IDOCalc.jl")
include("NCalc.jl")
using PyPlot
Al=logspace(-8,0,27)
w=logspace(-5,12,54)
yc=zeros(31,54)
ys=zeros(30,54)
Gc=zeros(31,54)
Gs=zeros(30,54)
cosacc=zeros(27,54)
sinacc=zeros(27,54)
cosstab=zeros(27,54)
sinstab=zeros(27,54)

for i=1:length(Al)
    (IDOC,IDOS)=IDOCalc(Al[i])
    for j=1:length(w)
        (yc[:,j],ys[:,j],Gc[:,j],Gs[:,j])=MP(1,w[j],FCMat,FCMatS,C_C,C_S,IDOC,IDOS)
        (cosacc[i,j],sinacc[i,j],cosstab[i,j],sinstab[i,j])=NCalc(1,w[j],FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm)
    end
    println(i)
end
