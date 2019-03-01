
using SymPy
using GenericSVD
using LinearAlgebra
function RUNC(NG,Al,n,d,cl,NoMo)
    # uses NG grid points (int)
    # uses Al for alpha - should be a big float
    # n for number of coarse grid pts on -1,1 (int)
    # d for degree of gram poly (int)
    # cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    # opt: 0 for even, 1 for odd
    # lam - weight of h1,h2
    # mu - weight of u
    
    


    include("ConvertToBig.jl")
    include("dftmtx.jl")
    include("Greens.jl")
    include("hsol.jl")
    include("coeffcalc.jl")
    include("Spaces.jl")
    include("All.jl")
    include("GreensInt.jl")
    include("FuncOut.jl")
    x=Sym("x")

    h1=2/(n-1);
    xcoarse1=collect(-1:h1:1-h1);

    h2=(2)/((n-1)*NG/(n-1));

    xx=collect(linspace(-1,1-h1,NG));
    
    f(x)=fmaker(d,n)
    ef=ConvertToBig(evalasarray(f(x),xx))
    #ef2=ConvertToBig(evalasarray((-1)^opt*f(x),xx))
    #ef=vcat(ef,ef2)
    #efff=ConvertToBig(evalasarray(f(x),xcoarse1))
    M=round(Int,round(Int,floor(NoMo/2)))
    println(M)
  
    ZZ=round(Int,((cl*NG/2+NG)))
    println(ZZ)
    LD=round(Int,NG/(n-1)+1)
    Modes=vcat(collect(1:M+1),collect(ZZ-M+1:ZZ))
    println(size(Modes))
    

    A=dftmtx(ZZ)
    A=adjoint(A)
    A=A/sqrt(ZZ)
    

    A1=A[1:NG,Modes]
    B=A[:,Modes]
    
    #ACont2=A[round(Int,ZZ/2)+1:round(Int,ZZ/2)+NG,:]

    
    LL=round(Int,ZZ/2)
    fd=2*((cl+2)/2);
    K=zeros(BigFloat,ZZ)
    K=-pi*1im*vcat(collect(0:LL-1),0,collect(-(LL-1):-1))

    
    # derivative coefficients pertaining to the FC
    qq=vcat(K[1:M+1],K[(ZZ)-M+1:(ZZ)])/fd;
    D = Diagonal(I,(length(qq)))-Al*Diagonal{BigFloat}(qq.*qq)
    D_Dag=pinv(D,1e-40)

    GG(x)=totalGreen(d,n,Al,xcoarse1[end])
    
    (h1,h2)=hsol(xx,Al,xcoarse1[end])


    EGGf=ConvertToBig(evalasarray(GG(x),xx))
    Eh1f=ConvertToBig(evalasarray(h1(x),xx))
    Eh2f=ConvertToBig(evalasarray(h2(x),xx))
    U_G=ConvertToBig(evalasarray(GG(x),xcoarse1))



    AugMat=A1;
   

 
    AugVec=ef[1:NG]


    (U,S,V)=GenericSVD.svd(AugMat)
    S=Diagonal{BigFloat}(S)
    S_Dag=pinv(S,1e-40)
    
    Res=V*(S_Dag*(transpose(U)*AugVec))

    
    
    

    Cin=Res
    
    F1=real(B*Cin);
    #FF1=real(ACont1*D_Dag*Cin)#+Res[end-1]*Eh1f+Res[end]*Eh2f);

acc=norm(F1[1:NG]-ef[1:NG])

    #FF2=real(ACont1*D_Dag*Cin+Res[end-1]*Eh1f+Res[end]*Eh2f)
    return F1,acc



    
  
end




