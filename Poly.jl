using SymPy
using GenericSVD
using LinearAlgebra
using PyPlot
function Poly(NG,Al,n,d,cl)
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

    (cobb,EGGf,Eh1f,Eh2f)=BCoeffCalc(xx,xcoarse1,Al,d,1)

    BB1=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f
    (cobb,EGGf,Eh1f,Eh2f)=BCoeffCalc(xx,xcoarse1,Al,d,2)
    BB2=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f
    (cobb,EGGf,Eh1f,Eh2f)=BCoeffCalc(xx,xcoarse1,Al,d,3)
    BB3=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f

    
    figure(1)

    plot(xcoarse,efff,label="Polynomial")
   
    # coarse grid things
    LD=round(Int,NG/(n-1))
  

    plot(xcoarse,BB1[1:LD:end],label="Func")
    plot(xcoarse,BB2[1:LD:end],"--",label="2nd Der")
    plot(xcoarse,BB3[1:LD:end],":",label="Sobolev")
    legend()
    

end
