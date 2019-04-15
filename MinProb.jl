using SymPy
using LinearAlgebra
using PyPlot
using GenericSVD
function MinProb(d,Al,Lam,Mu,Gam)


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

    x=symbols("x")


    setprecision(230);
    D=BigFloat(350)/BigFloat(100);
    A=BigFloat(125)/BigFloat(100);
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
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)

    GG(x)=TotalPoly(d,n,Al);
    EGGf=evalasarray(GG(x),xx);
    (h1,h2)=hsol(xx,Al,0)
    Eh1f=(evalasarray(h1(x),xx))
    Eh2f=(evalasarray(h2(x),xx))

    if d==0
    c1=-0.000014098152457
    elseif d==1
    c1=0.000009248317084
    elseif d==2
    c1=0.001398898944611
    elseif d==3
    c1=-0.018604242523930
    elseif d==4
    c1=-0.126012962108665
    elseif d==5
    c1=-0.555408319762897
    elseif d==6
    -1.958059513322267
    end


    BF=EGGf-(-1)^d*c1*Eh1f-c1*Eh2f;

    (M,B,C)=CosMtx(D,A,n,F);

    (s1,s2)=size(B);


    (D1,D2)=CosDer2(D,A,n,F);

    DO=DiffOp(D,A,n,F,Al,1);
    IDO=zeros(BigFloat,size(DO))
    for i=1:length(DO[1,:])
        IDO[i,i]=1/DO[i,i];
    end
    Mat1=B*IDO;
    r=setdiff(1:s1,1:10:s1);
    Mat1=Mat1[setdiff(1:end,r),:]

    Ent_1=hcat(B,zeros(BigFloat,s1,2));
    Ent_2=hcat(Mu*Mat1,Mu*Al*Eh1f[1:10:end],Mu*Al*Eh2f[1:10:end]);
    Ent_3=hcat(zeros(BigFloat,s1,s2),Lam*Al*Eh1f,Lam*Al*Eh2f)
    Ent_4=hcat(Gam*B,zeros(BigFloat,s1,2));

    MM=vcat(Ent_1,Ent_2,Ent_3,Ent_4);
    RHS=vcat(y,BF[1:10:end],0*Eh1f,0*y)
    #(Um,Sm,Vm)=svd(MM)
    #SmD=pinv(Sm,1e-40);
    #SmD=Diagonal{BigFloat}(SmD')
    #Res=Vm*(SmD*(Um'*(RHS)));
    Res=MM\RHS;


    return M*Res[1:end-2],Res[end-1:end]








    end