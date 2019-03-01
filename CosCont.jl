using SymPy
using GenericSVD
using LinearAlgebra

include("CosMtx.jl")
include("CosDer.jl")
include("All.jl")
include("ConvertToBig.jl")
include("hsol.jl")
include("Spaces.jl")
include("dftmtx.jl")


x=Sym("x")

setprecision(100);

d=0;
D=3.5;
A=1.25;
n=10;
F=10;

f(x)=fmaker(d,n+1);
xx=linspace(-1,1,F*(n-1)+1);
y=ConvertToBig(evalasarray(f(x),xx));


x_g=collect(A:1/(n-1):D-A);

(M,B)=CosMtx(D,A,n,F);
(D1,D2)=CosDer(D,A,n,F);

(U,S,V)=GenericSVD.svd(B);
S=Diagonal(S);
Sd=pinv(S,1e-50);
Res1=V*(Sd*(U'*y));
plot(A*Res1)




