using SymPy
using LinearAlgebra
using PyPlot
function SinCont(d)


    include("SinMtx.jl")
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
    xx=linspace(BigFloat(-1),BigFloat(1),Int(F*(n-1)+1));
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)

    (M,B,C)=SinMtx(D,A,n,F);

    (s1,s2)=size(B);


    (Uc,Sc,Vc)=GenericSVD.svd(B);
    Sc=Diagonal{BigFloat}(Sc);
    Scd=pinv(Sc,1e-40);
    Res1=Vc*(Scd*(Uc'*y));

    fc=M*Res1;
    f_AD=fc[428:428+90];
    uc=M*(IDO*Res1[1:end-2])
    u_AD=uc[428:428+90]

    return maximum(abs.(f_AD-y)),Res1,fc

end
