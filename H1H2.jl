using SymPy
using LinearAlgebra
using PyPlot
using DelimitedFiles


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
    d=0
    Al=.01

x=Sym("x")

setprecision(200);

D=3.5;
A=1.25;
n=10;
    F=10;

    hc=1/(n-1);
    Coarse=collect(BigFloat,-D:hc:D-hc)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(BigFloat,-D:hf:D-hf)
    F1=collect(BigFloat,A:hf:D-A);

h=1/(n-1);
h1=h/F;


f(x)=fmaker(d,n+1);
    xx=linspace(-1,1-2/n,F*(n-1)+1);
    y=zeros(BigFloat,size(xx))
    y=evalasarray(f(x),xx)
    figure(10)


x_g=collect(BigFloat,A:1/(n-1)/F:D-A);
    xf=collect(BigFloat,-D:h1:D-h1);

    GG(x)=totalGreen(d,n,Al,1-2/n);
    EGGf=evalasarray(GG(x),xx);
    (h1,h2)=hsol(xx,Al,1-2/n)
    Eh1f=(evalasarray(h1(x),xx))
    Eh2f=(evalasarray(h2(x),xx))
plot(xx,EGGf,xx,Eh1f,xx,Eh2f)

M=hcat(Eh1f,Eh2f,EGGf)

open("H1H2G.txt","w") do f
    writedlm(f,M)
end

    


