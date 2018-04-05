using SymPy

include("All.jl")
include("hsol.jl")
include("dftmtx.jl")
include("GreensInt.jl")
include("Greens.jl")
include("ConvertToBig.jl")
x=Sym("x")
Al=.1
NG=32
n=9
cl=2
d=2

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
EGGf=ConvertToBig(evalasarray(GG(x),xx))
Eh1f=ConvertToBig(evalasarray(h1(x),xx))
Eh2f=ConvertToBig(evalasarray(h2(x),xx))



Z=[(Eh1'*Eh1) (1im*EH1'*Eh1);(Eh2'*EH2) (1im*EH2'*EH2)]
zz=[0;0]
co=Z\zz

println(co)

