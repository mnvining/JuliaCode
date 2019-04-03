using SymPy
using LinearAlgebra
using PyPlot
function CosContLM(d,Al,Lam,Mu,N)
    # N is number of total points, divisible by 9!!


    include("CosMtx.jl")
    include("CosDer.jl")
    include("CosDer2.jl")
    include("All.jl")
    include("ConvertToBig.jl")
    include("hsol.jl")
    include("Spaces.jl")
    include("dftmtx.jl")
    include("coeffcalc.jl")
    include("DiffOp.jl")
    include("GreensInt.jl")
    include("myhouse.jl")

    x=Sym("x")


    setprecision(200);
    CL=BigFloat(N)/BigFloat(10);
    D=BigFloat(CL)/BigFloat(2);
    A=BigFloat(1)/BigFloat(2)*(D-1);
    n=10;
    F=10;

    hc=BigFloat(1)/BigFloat(n-1);
    Coarse=collect(BigFloat,-D+hc:hc:D)
    LC=length(Coarse)
    hf=BigFloat(1)/BigFloat((n-1)*F);
    Fine=collect(BigFloat,-D+hf:hf:D)
    F1=collect(BigFloat,A:hf:D-A);


    f(x)=fmaker(d,n+1);
    xx=linspace(BigFloat(-1),BigFloat(0),Int(F*(n-1)+1));
    x2=linspace(BigFloat(-1),BigFloat(0),Int(n));
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)

    GG(x)=totalGreen(d,n+1,Al,0);
    EGGf=evalasarray(GG(x),xx);
    EG2=evalasarray(GG(x),x2);
    println(EGGf[1:F:end]-EG2)
    G2=evalasarray(GG(x),xx);
    G2[1:F:end]=zeros(n,1);
    GC=EGGf-G2;
    (h1,h2)=hsol(xx,Al,0)
    Eh1f=(evalasarray(h1(x),xx))
    T1=evalasarray(h1(x),xx)
    T1[1:F:end]=zeros(n,1);
    E1c=Eh1f-T1;
    Eh2f=(evalasarray(h2(x),xx))
    T2=evalasarray(h2(x),xx)
    T2[1:F:end]=zeros(n,1);
    E2c=Eh2f-T2;

    (M,B,C)=CosMtx(D,A,n,F);

    (s1,s2)=size(B);

    DO=DiffOp(D,A,n,F,Al,1);
    IDO=zeros(BigFloat,size(DO))
    for i=1:length(DO[1,:])
        IDO[i,i]=1/DO[i,i];
    end





    Ent_1=hcat(B,zeros(BigFloat,s1,2));
    #Ent_2=hcat(B*IDO,Eh1f,Eh2f);
    Ent_3=hcat(zeros(BigFloat,s1,s2),Lam*Eh1f,Lam*Eh2f)
    P=B*IDO
    P2=B*IDO
    P2[1:F:end,:]=zeros(BigFloat,length(1:F:s1),s2)
    US=P-P2;# this is the coarse grid eval of u I THINK
    Ent_4=hcat(US,E1c,E2c)


    #AugMat=vcat(Ent_1,Ent_2,Ent_3,Ent_4)
    AugMat=vcat(Ent_1,Ent_4,Ent_3)
    #AugVec=vcat(y,EGGf,zeros(BigFloat,size(Eh1f)),zeros(BigFloat,size(Eh1f)))
    AugVec=vcat(y,GC,zeros(BigFloat,size(Eh1f)))

    (Uc,Sc,Vc)=GenericSVD.svd(AugMat);
    Sc=Diagonal{BigFloat}(Sc);
    Scd=pinv(Sc,1e-40);
    Res1=Vc*(Scd*(Uc'*AugVec));
    LL=Int(428)
    fc=M*Res1[1:end-2];
    f_AD=fc[LL:LL+90];
    uc=M*(IDO*Res1[1:end-2])
    u_AD=uc[LL:LL+90]

    return norm(f_AD-y),norm(u_AD[1:10:end])/norm(y[1:10:end]),fc,uc,Res1[end-1:end]

end
