using GenericSVD

function MP(d,w,FCMat,FCMatS,C_C,C_S,IDOC,IDOS)
ffc=copy(FCMat[:,d+1])
ffs=copy(FCMatS[:,d+1])

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)

MatC=vcat(IDOC'*(IDOC*Vc),w*Sc_Full);
MatS=vcat(IDOS'*(IDOS*Vs),w*Ss_Full);

Rc=vcat(IDOC*ffc,zeros(BigFloat,size(IDOC*ffc)));
Rs=vcat(IDOS*ffs,zeros(BigFloat,size(IDOS*ffs)));

Gc=MatC\Rc;
Gs=MatS\Rs;

yc=C_C*IDOC*(ffc-Vc*Gc);
ys=C_S*IDOS*(ffs-Vs*Gs;
return yc,ys
end
