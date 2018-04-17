using PyPlot
using SymPy
using GenericSVD
function RUN2(NG,Al,n,d,cl,opt,plotopt,c_1,c_2)
    # uses NG grid points (int)
    # uses Al for alpha - should be a big float
    # n for number of coarse grid pts on -1,1 (int)
    # d for degree of gram poly (int)
    # cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    # opt: 1 for func, 2 for 2nd der, 3, for sobolev
    # c_1 - coeff of f c_2 coeff d^2 f
    
    


include("ConvertToBig.jl")
include("dftmtx.jl")
include("Greens.jl")
include("hsol.jl")
include("coeffcalc.jl")
include("All.jl")
include("GreensInt.jl")
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
    #println(length(Modes))

A=createFWDcont(NG,Modes,cl)
B=createBackDoubleDomain(NG,Modes,cl)

GG(x)=totalGreen(d,n,Al)
d2GG(x)=diff(GG(x),x,2)
(h1,h2)=hsol(xx,Al)
d2h1(x)=diff(h1(x),x,2)
    d2h2(x)=diff(h2(x),x,2)
    GPrime(x)=diff(GG(x),x,1)
    h1Prime(x)=diff(h1(x),x,1)
    h2Prime(x)=diff(h2(x),x,1)
    EGG=ConvertToBig(evalasarray(GG(x),xcoarse1))
    Eh1=ConvertToBig(evalasarray(h1(x),xcoarse1))
    Eh2=ConvertToBig(evalasarray(h2(x),xcoarse1))
    eh1pr=ConvertToBig(evalasarray(h1Prime(x),xcoarse1))
    eh2pr=ConvertToBig(evalasarray(h2Prime(x),xcoarse1))
    Gpr=ConvertToBig(evalasarray(GPrime(x),xcoarse1))
    EG2=ConvertToBig(evalasarray(d2GG(x),xcoarse1))
    EH1=ConvertToBig(evalasarray(d2h1(x),xcoarse1))
    EH2=ConvertToBig(evalasarray(d2h2(x),xcoarse1))
    if opt==1
        Z=[Eh1'*Eh1 Eh2'*Eh1; Eh1'*Eh2 Eh2'*Eh2]
zz=[-EGG'*Eh1; -EGG'*Eh2]
    elseif opt==2
        Z=[EH1'*EH1 EH2'*EH1; EH1'*EH2 EH2'*EH2]
zz=[-EG2'*EH1; -EG2'*EH2]
    elseif opt==3
        Z=[(Eh1'*Eh1+eh1pr'*eh1pr) (Eh1'*Eh2+eh1pr'*eh2pr); (Eh1'*Eh2+eh1pr'*eh2pr) (Eh2'*Eh2+eh2pr'*eh2pr)]
        zz=[-EGG'*Eh1-Gpr'*eh1pr; -EGG'*Eh2-Gpr'*eh2pr]
    elseif opt==4
        Z=[(Eh1'*Eh1+eh1pr'*eh1pr+EH1'*EH1) (Eh1'*Eh2+eh1pr'*eh2pr+EH1'*EH2); (Eh1'*Eh2+eh1pr'*eh2pr+EH1'*EH2) (Eh2'*Eh2+eh2pr'*eh2pr+EH2'*EH2)]
        zz=[-EGG'*Eh1-Gpr'*eh1pr-EG2'*EH1; -EGG'*Eh2-Gpr'*eh2pr-EG2'*EH1]
    else
        

        Z=[(c_1^2*Eh1'*Eh1+c_1*c_2*Eh1'*EH1+c_1*c_2*EH1'*Eh1+c_1*c_2*Eh1'*EH1+c_2^2*EH1'*EH1) (c_1^2*Eh2'*Eh1+c_1*c_2*Eh2'*EH1+c_1*c_2*EH2'*EH1+c_1*c_2*Eh2'*EH1+c_2^2*EH2'*EH1);(c_1^2*Eh1'*Eh2+c_1*c_2*Eh1'*EH2+c_1*c_2*EH1'*Eh2+c_1*c_2*Eh1'*EH2+c_2^2*EH1'*EH2) (c_1^2*Eh2'*Eh2+c_1*c_2*Eh2'*EH2+c_1*c_2*EH2'*EH2+c_1*c_2*Eh2'*EH2+c_2^2*EH2'*EH2)]
        zz=[-c_1^2*EGG'Eh1-c_1*c_2*EGG'*EH1-c_1*c_2*EG2'*Eh1-c_1*c_2*EGG'*Eh1-c_2^2*EG2'*EH1;-c_1^2*EGG'Eh2-c_1*c_2*EGG'*EH2-c_1*c_2*EG2'*Eh2-c_1*c_2*EGG'*Eh2-c_2^2*EG2'*EH2]
        end

    EGGf=ConvertToBig(evalasarray(GG(x),xx))
    Eh1f=ConvertToBig(evalasarray(h1(x),xx))
    Eh2f=ConvertToBig(evalasarray(h2(x),xx))

    cobb=Z\zz


    BB=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f

D=createId2Derivative(NG,M,Al,cl)
D_Dag=pinv(D,1e-14)

AugMat=vcat(A,A*D_Dag')

AugVec=vcat(ef,BB)



Cin=AugMat\AugVec


    FF1ex=real(B*D_Dag*Cin)
    FF1=real(A*D_Dag*Cin)
    FF2=real(A*Cin)
    FF2ex=real(B*Cin)
acc=norm(FF1-BB)
    acc2=norm(FF2-ef)
    # coarse grid things
    LD=round(Int,NG/(n-1))
    solcoarse=FF1[1:LD:end];

    
    stab=norm(solcoarse)/norm(efff)
    #figure(d+1)
    if plotopt==1
        figure()

        plot(xx,FF1,xx,FF2,xx,ef,xx,BB,xx,EGGf)
        legend(("FullOp","Fc","func","Bf","Greens"))
    elseif plotopt==2
        plot(xx,ef,xx,BB,xx,EGGf)
        legend(("Func","Bf","Greens"))
    end
    
    

    
    NB=round(Int,floor(length(BB)/2))

    return acc,acc2,stab,BB[NB]
end




