# Run only after running Storage.jl
using GenericSVD
fc=FCMat[:,1]
fs=FCMatS[:,1]

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)
w=1
tol=1e-16 # working tolerance for sing value contribution

Pc=Vc'*(IDOC.^2*(Vc))+w*(1/tol^2)*Sc_Full.^2;
Ps=Vs'*(IDOS.^2*(Vs))+w*(1/tol^2)*Ss_Full.^2;
uc=-IDOC*fc;
us=-IDOS*fs;
qc=-Vc'*IDOC'*uc;
qs=-Vs'*IDOS'*us;

# Solve p*x=-q for KKT conditions
Gc=Pc\(-qc)
Gs=Ps\(-qs)
