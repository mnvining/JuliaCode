using GenericSVD

function MP(d,w,FCMatC,FCMatS,C_C,C_S,IDOC,IDOS,tol)
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

    

    MatC=vcat((C_C[1:10:end,:]*IDOC*Vc),w*Sc_Full);
    MatS=vcat((C_S[1:10:end,:]*IDOS*Vs),w*Ss_Full);

    Rc=vcat(C_C[1:10:end,:]*IDOC*ffc,zeros(BigFloat,size(IDOC*ffc)));
    Rs=vcat(C_S[1:10:end,:]*IDOS*ffs,zeros(BigFloat,size(IDOS*ffs)));

    if tol == 0
        Gc=MatC\Rc;
        Gs=MatS\Rs;
    else
        (U,S,V)=GenericSVD.svd(MatC);
        (U2,S2,V2)=GenericSVD.svd(MatS);
        S=Diagonal{BigFloat}(S);
        Sd=pinv(S,tol);
        S2=Diagonal{BigFloat}(S2);
        S2d=pinv(S2,tol);

        Gc=V*Sd*U'*Rc;
        Gs=V2*S2d*U2'*Rs;
    end



    yc=(ffc-Vc*Gc);
    ys=(ffs-Vs*Gs);
    return yc,ys,Gc,Gs
end
