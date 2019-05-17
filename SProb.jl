# Run only after running Storage.jl
using GenericSVD
function SProb(d,tol)
fc=FCMat[:,d+1]
fs=FCMatS[:,d+1]

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)
w=0
#
Pc=Vc'*(IDOC.^2*(Vc))+w*(1/tol^2)*(Sc_Full'*Sc_Full);
Ps=Vs'*(IDOS.^2*(Vs))+w*(1/tol^2)*(Ss_Full'*Ss_Full);
uc=IDOC*fc;
us=IDOS*fs;
qc=-Vc'*IDOC'*uc;
qs=-Vs'*IDOS'*us;

# Solve p*x=-q for KKT conditions
Gc=Pc\(-qc)
Gs=Ps\(-qs)

# reconstruct with a -!
yc=C_C*IDOC*fc-C_C*IDOC*Gc;
ys=C_S*IDOS*fs-C_S*IDOS*Gs;
return yc,ys
end
