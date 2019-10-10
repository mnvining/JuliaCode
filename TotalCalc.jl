using GenericSVD
function TotalCalc(d,w,tol,FCMatC,FCMatS,C_C,C_S,IDOC,IDOS,m,mm)
    ffc=copy(FCMatC[:,d+1])
    ffs=copy(FCMatS[:,d+1])

    setprecision(230);
    D=BigFloat(350)/BigFloat(100);
    A=BigFloat(125)/BigFloat(100);
    n=10;
    F=10;
    x=Sym("x")


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

    #Gc=MatC\Rc;
    #Gs=MatS\Rs;

    (U,S,V)=GenericSVD.svd(MatC);
    (U2,S2,V2)=GenericSVD.svd(MatS);
    S=Diagonal{BigFloat}(S);
    Sd=pinv(S,tol);
    S2=Diagonal{BigFloat}(S2);
    S2d=pinv(S2,tol);

    Gc=V*Sd*U'*Rc;
    Gs=V2*S2d*U2'*Rs;



    yc=(ffc-Vc*Gc);
    ys=(ffs-Vs*Gs);

    f(x)=fmaker(d,n);
    xx=linspace(BigFloat(-1),BigFloat(1),Int(F*(n-1)+1));
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)

    CAcc=maximum(abs.(C_C*yc.-y));
    SAcc=maximum(abs.(C_S*ys.-y));

    uc=C_C*IDOC*yc;
    us=C_S*IDOS*ys;
    CStab=norm(uc[1:10:end])/norm(y[1:10:end])
    SStab=norm(us[1:10:end])/norm(y[1:10:end])

    return CStab,SStab,CAcc,SAcc, maximum(abs.(m*yc)),maximum(abs.(mm*ys))
end
