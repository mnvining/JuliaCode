using SymPy
using LinearAlgebra
using PyPlot
function CosCont(d)


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

    (M,B,C)=CosMtx(D,A,n,F);

    (s1,s2)=size(B);


    (Uc,Sc,Vc)=GenericSVD.svd(B);
    Sc=Diagonal{BigFloat}(Sc);
    Scd=pinv(Sc,1e-40);
    Res1=Vc*(Scd*(Uc'*y));

    fc=M*Res1[1:end-2];
    f_AD=fc[428:428+90];
    uc=M*(IDO*Res1[1:end-2])
    u_AD=uc[428:428+90]

    return norm(f_AD-y),norm(u_AD[1:10:end])/norm(y[1:10:end]),Res1

end
