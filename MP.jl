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



MatC=vcat((IDOC*Vc),w*Sc_Full);
MatS=vcat((IDOS*Vs),w*Ss_Full);
println(typeof(MatC))

Rc=vcat(IDOC*ffc,zeros(BigFloat,size(IDOC*ffc)));
Rs=vcat(IDOS*ffs,zeros(BigFloat,size(IDOS*ffs)));

Gc=MatC\Rc;
Gs=MatS\Rs

println(rank(IDOC))
println(size(IDOC))

yc=(ffc-Vc*Gc);
ys=(ffs-Vs*Gs);
return yc,ys
end
