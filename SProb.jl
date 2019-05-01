# Run only after running Storage.jl
using GenericSVD

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)
w=1
tol=1e-16 # working tolerance for sing value contribution
P=V'*(D_Dag.^2*(V))+w*(1/tol^2)*S_Full.^2;
u=-D_Dag*f_c;
q=-V'*D_Dag'*u;

# Solve p*x=-q for KKT conditions

G=P\q
