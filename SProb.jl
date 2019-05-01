# Run only after running Storage.jl
using GenericSVD

(U,S,V)=GenericSVD.svd(C)
S_Full=Diagonal{BigFloat}(S)
w=1
tol=1e-16 # working tolerance for sing value contribution
P=V'*(D_Dag.^2*(V))+w*(1/tol^2)*S_Full.^2;
u=-D_Dag*f_c;
q=-V'*D_Dag'*u;

# Solve p*x=-q for KKT conditions

G=P\q
