using SymPy
using LinearAlgebra
using PyPlot
function MinProb(d,Al,Lam)


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

    GG(x)=totalGreen(d,n+1,Al,0);
    EGGf=evalasarray(GG(x),xx);
    (h1,h2)=hsol(xx,Al,0)
    Eh1f=(evalasarray(h1(x),xx))
    Eh2f=(evalasarray(h2(x),xx))

    c1=CCalc(d,Al)

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
    Mat1=Mat1[setdiff(1:end,r),:];
    Z1=zeros(BigFloat,10,10);
    P=B'*B-Lam*Diagonal{BigFloat}(I,s2);
    println(size(P))
    M1=hcat(P,Mat1')
    M2=hcat(Mat1,Z1)
    MM=vcat(M1,M2)








    end
