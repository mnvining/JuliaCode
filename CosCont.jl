using SymPy
using LinearAlgebra
using PyPlot
function CosCont(d,Al,Lam)


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


    #x_g=collect(BigFloat,A:1/(n-1)/F:D-A);
    #xf=collect(BigFloat,-D-h1:h1:D);

    GG(x)=totalGreen(d,n+1,Al,0);
    EGGf=evalasarray(GG(x),xx);
    (h1,h2)=hsol(xx,Al,0)
    Eh1f=(evalasarray(h1(x),xx))
    Eh2f=(evalasarray(h2(x),xx))
   

    (M,B,C)=CosMtx(D,A,n,F);

    
    println(size(B))

    
    (s1,s2)=size(B);

    
    (D1,D2)=CosDer2(D,A,n,F);

    DO=DiffOp(D,A,n,F,Al,1);
    IDO=zeros(BigFloat,size(DO))
    for i=1:length(DO[1,:])
        IDO[i,i]=1/DO[i,i];
    end
    
    
    Ent_1=hcat(B,zeros(BigFloat,s1,2));
    Ent_2=hcat(B*IDO,Eh1f,Eh2f);
    Ent_3=hcat(zeros(BigFloat,s1,s2),Lam*Eh1f,Lam*Eh2f)
    
    AugMat=vcat(Ent_1,Ent_2,Ent_3)
    AugVec=vcat(y,EGGf,zeros(BigFloat,size(Eh1f)))
    #figure(88)
    #plot(y)

    #(Res1,PP)=SolveViaQR(AugMat,AugVec)
    #(Res1,PP)=SolveViaQR(B,y);
    (Uc,Sc,Vc)=GenericSVD.svd(AugMat);
    Sc=Diagonal{BigFloat}(Sc);
    Scd=pinv(Sc,1e-40);
    Res1=Vc*(Scd*(Uc'*AugVec));
    #Res1=AugMat\AugVec
    

    
    #Res1=B\y;

    fc=M*Res1[1:end-2];
    f_AD=fc[428:428+90];
    uc=M*(IDO*Res1[1:end-2])
    u_AD=uc[428:428+90];    
    
    

    return norm(f_AD-y),uc,fc,y#,AugMat,AugVec

end

