using SymPy
using LinearAlgebra
using PyPlot
using GenericSVD
function CosContNoD(d,Al)


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

    (M,B,C)=CosMtx(D,A,n,F);

    (U,S,V)=GenericSVD.svd(B);
    S=Diagonal{BigFloat}(S);
    SD=GenericSVD.pinv(S,1e-40);
    Res1=V*(SD*(U'*y))
    


    fc=M*Res1;
    f_AD=fc[428:428+90];

    DO=DiffOp(D,A,n,F,Al,1);
    IDO=zeros(BigFloat,size(DO))
    for i=1:length(DO[1,:])
        IDO[i,i]=1/DO[i,i];
    end

    uc=M*IDO*Res1;
    u_AD=uc[428:428+90]
    
    
    
    

    return norm(f_AD-y),fc,y,uc#,AugMat,AugVec

end

