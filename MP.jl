using GenericSVD

function MP(d,w,FCMatC,FCMatS,C_C,C_S,IDOC,IDOS)
ffc=copy(FCMatC[:,d+1])
ffs=copy(FCMatS[:,d+1])

(Uc,Sc,Vc)=GenericSVD.svd(C_C)
(Us,Ss,Vs)=GenericSVD.svd(C_S)
Sc_Full=Diagonal{BigFloat}(Sc)
Ss_Full=Diagonal{BigFloat}(Ss)
(s1,s2)=size(Sc_Full);
(s11,s12)=size(Ss_Full)
Sc_Full=copy(Sc_Full)+zeros(BigFloat,s1,s2)
Ss_Full=copy(Ss_Full)+zeros(BigFloat,s11,s12)



MatC=vcat((C_C[1:10:end]*IDOC*Vc),w*Sc_Full);
MatS=vcat((C_S[1:10:end]*IDOS*Vs),w*Ss_Full);

Rc=vcat(C_C[1:10:end]*IDOC*ffc,zeros(BigFloat,size(IDOC*ffc)));
Rs=vcat(C_S[1:10:end]*IDOS*ffs,zeros(BigFloat,size(IDOS*ffs)));

Gc=MatC\Rc;
Gs=MatS\Rs

yc=(ffc-Vc*Gc);
ys=(ffs-Vs*Gs);
return yc,ys
end
