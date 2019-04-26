using SymPy
using LinearAlgebra
using PyPlot
using GenericSVD
using Roots
function ZeroTest(d,Al,Lam,Mu,Gam)


    include("CosMtx2.jl")
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
    CC=symbols("CC")


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


    f(x)=fmaker(d,n);
    println(f(x))
    xx=linspace(BigFloat(-1),BigFloat(1),Int(F*(n-1)+1));
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)
    ya=vcat(1,0*y,1)

    #if d==0
    #c1=-1.409815245715622e-05
    #elseif d==1
    #c1=9.248317083713627e-06
    #elseif d==2
    #c1=0.001398898944611
    #elseif d==3
    #c1=-0.018604242523930
    #elseif d==4
    #c1=-0.126012962108665
    #elseif d==5
    #c1=-0.555408319762897
    #elseif d==6
    #-1.958059513322267
    #end

    (M,B,C)=CosMtx2(D,A,n,F);

    (U,S,V)=svd(B);
    S=Diagonal{BigFloat}(S)
    Sdag=pinv(S,1e-40);
    Res=V*(Sdag*(U'*ya));

    CF=M*Res;

    (D1,D2)=CosDer2(D,A,n,F);

    DO=DiffOp(D,A,n,F,Al,1);
    IDO=zeros(BigFloat,size(DO))
    for i=1:length(DO[1,:])
        IDO[i,i]=1/DO[i,i];
    end
    CU=M*IDO*Res;
    UG=CU[428:518]

    FG=CF[428:518]


    return CF,norm(FG.-0*y),norm(UG[1:10:end])/norm(y[1:10:end]),Coarse








    end
