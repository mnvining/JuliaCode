# Run only after running Storage.jl
using GenericSVD
function SProb(d,tol,FCMat,FCMatS,C_C,C_S,IDOC,IDOS)
ffc=copy(FCMat[:,d+1])
ffs=copy(FCMatS[:,d+1])

#plot(C_S*IDOS*fs)

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)
w=1

Pc=Vc'*(IDOC.^2*(Vc))+w*(1/tol^2)*(Sc_Full'*Sc_Full);
Ps=Vs'*(IDOS.^2*(Vs))+w*(1/tol^2)*(Ss_Full'*Ss_Full);
uc=IDOC*ffc;
us=IDOS*ffs;
qc=Vc'*IDOC'*uc;
qs=Vs'*IDOS'*us;

# Solve p*x=-q for KKT conditions
Gc=Pc\(qc)
Gs=Ps\(qs)

semilogy(2*Vc'*ffc)

# reconstruct with a -!
yc=C_C*IDOC*ffc-C_C*IDOC*Gc;
ys=C_S*IDOS*ffs-C_S*IDOS*Gs;
return yc,ys
end
