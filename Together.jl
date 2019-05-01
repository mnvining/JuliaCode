# Run only after running Storage.jl
using GenericSVD
fc=FCMat[:,1]
fs=FCMatS[:,1]
f_c=vcat(fc,fs)
C=hcat(C_C,C_S)

Ent_1=hcat(IDOC,zeros(31,30))
Ent_2=hcat(zeros(30,31),IDOS)
D_Dag=vcat(Ent_1,Ent_2)

(U,S,V)=GenericSVD.svd(C)
S_Full=Diagonal{BigFloat}(S);

w=1
tol=1 # working tolerance for sing value contribution

P=V'*(D_Dag.^2*(V))+w*(1/tol^2)*S_Full.^2;
u=-D_Dag*f_c;

q=-V'*D_Dag'*u;

# Solve p*x=-q for KKT conditions
G=P\(-q)
