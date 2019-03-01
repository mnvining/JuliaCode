using SymPy
using GenericSVD
using LinearAlgebra
using PyPlot
function RUNGraph(NG,Al,n,d,cl,opt)
    # uses NG grid points (int)
    # uses Al for alpha - should be a big float
    # n for number of coarse grid pts on -1,1 (int)
    # d for degree of gram poly (int)
    # cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    # opt: 1 for func, 2 for 2nd der, 3, for sobolev
    include("ConvertToBig.jl")
    include("dftmtx.jl")
    include("Greens.jl")
    include("hsol.jl")
    include("coeffcalc.jl")
    include("All.jl")
    include("GreensInt.jl")
    include("Spaces.jl")
    include("FuncOut.jl")
    include("BCoeffs.jl")
    
    x=Sym("x")

    xx=linspace(-1,1,NG+1)
    xx=xx[1:end-1]
    h=xx[2]-xx[1]
    x2=collect(-1:h:1+cl-h)
    xcoarse1=collect(linspace(-1,1,n))
    xcoarse=xcoarse1[1:end-1]
    hh=xcoarse[2]-xcoarse[1]
    xcoarseext=collect(-1:hh:1+cl-hh)
    f(x)=fmaker(d,n)
    ef=ConvertToBig(evalasarray(f(x),xx))
    efff=ConvertToBig(evalasarray(f(x),xcoarse))
    M=round(Int,floor(length(xcoarseext)/2))
  
    ZZ=Int(cl*NG/2+NG);
    Modes=vcat(collect(1:M+1),collect(ZZ-M+1:ZZ))

    A=createFWDcont(NG,Modes,cl)
    B=createBackDoubleDomain(NG,Modes,cl)

    (cobb,EGGf,Eh1f,Eh2f)=BCoeffCalc(xx,xcoarse1,Al,d,opt)

    BB=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f

    D=createId2Derivative(NG,M,Al,cl)
    D_Dag=pinv(D,1e-14)

    AugMat=vcat(A,A*D_Dag')

    AugVec=vcat(ef,BB)

    Cin=AugMat\AugVec

    Uex=real(B*D_Dag*Cin)
    U1=real(A*D_Dag*Cin)
    U2=real(A*Cin)
    U2ex=real(B*Cin)
    acc=norm(U1-BB)
    acc2=norm(U2-ef)
    # coarse grid things
    LD=round(Int,NG/(n-1))
    solcoarse=U1[1:LD:end];
  
    stab=norm(solcoarse)/norm(efff)
    stab2=norm(BB[1:LD:end])/norm(efff)

    NB=round(Int,floor(length(BB)/2))

    return stab,cobb[1],cobb[2],acc,stab2
end
