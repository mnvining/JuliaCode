# This works for what we want it to do 1/29/18
# Start coding in BigFloat precision to see if we get better results


using PyPlot
using SymPy
using GenericSVD
function RUN(NG,Al,n,d,cl)
# uses NG grid points
# uses Al for alpha
# n for number of coarse grid pts on -1,1
# d for degree of gram poly
# cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
# hard coded for debugging



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
f(x)=fmaker(d,n)
ef=evalasarray(f(x),xx)
M=round(Int,NG/4)
ZZ=Int(cl*NG/2+NG);
Modes=vcat(collect(1:M+1),collect(ZZ-M+1:ZZ))

A=createFWDcont(NG,Modes,cl)
B=createBackDoubleDomain(NG,Modes,cl)
# here's where we a) modify G
# b) modify the matrix we're doing the solve for in order to get the [C;CD*]\[f;Bf]
GG(x)=totalGreen(d,n,Al)
d2GG(x)=diff(GG(x),x,2)
(h1,h2)=hsol(xx,Al)
d2h1(x)=diff(h1(x),x,2)
d2h2(x)=diff(h2(x),x,2)
EGG=evalasarray(GG(x),xx)
Eh1=evalasarray(h1(x),xx)
Eh2=evalasarray(h2(x),xx)
EG2=evalasarray(d2GG(x),xx)
EH1=evalasarray(d2h1(x),xx)
EH2=evalasarray(d2h2(x),xx)

Z=[EH1'*EH1 EH2'*EH1; EH1'*EH2 EH2'*EH2]
zz=[-EG2'*EH1; -EG2'*EH2]

cobb=Z\zz
#println(cobb)

BB=EGG+cobb[1]*Eh1+cobb[2]*Eh2
println(norm(BB))

D=createId2Derivative(NG,M,Al,cl)
D_Dag=pinv(D,1e-14)

AugMat=vcat(A,A*D_Dag')

AugVec=vcat(ef,BB)



Cin=AugMat\AugVec


FF1=real(A*D_Dag*Cin)
FF2=real(A*Cin)
acc=norm(FF1-BB)
acc2=norm(FF2-ef)



    plot(xx,FF1,xx,BB)

    return acc,acc2

end


#title("Fourier Continuation of f1(x) with N=32, Modes=17")
