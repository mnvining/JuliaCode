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
    D=BigFloat(4);
    A=BigFloat(52)/BigFloat(35);
    n=10;
    F=10;

    hc=BigFloat(4)/BigFloat(35);
    Coarse=collect(BigFloat,-D+hc:hc:D)
    #Coarse=Coarse.-BigFloat(5)/BigFloat(90)
    #LC=length(Coarse)
    hf=BigFloat(hc)/BigFloat(F);
    Fine=collect(BigFloat,-D+hc:hf:D)
    #Fine=Fine.-(BigFloat(5)/BigFloat(900));
    F1=collect(BigFloat,A:hf:D-A);



    f(x)=fmaker(d,n);
    xx=linspace(BigFloat(-1),BigFloat(1),Int(F*(n-1)+1));
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)
    if d==9
        y=y/3
    end

    (M,B,C)=SinMtx(D,A,n,F);

    (s1,s2)=size(B);


    (Uc,Sc,Vc)=GenericSVD.svd(B);
    Sc=Diagonal{BigFloat}(Sc);
    Scd=pinv(Sc,1e-40);
    Res1=Vc*(Scd*(Uc'*y));

    fc=M*Res1;
    f_AD=fc[428:428+90];


    return maximum(abs.(f_AD-y)),Res1,fc

end
