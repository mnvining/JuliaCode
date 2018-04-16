using PyPlot
using SymPy
using GenericSVD
include("TryinBig.jl")
include("All.jl")
include("dftmtx.jl")
include("coeffcalc.jl")
include("hsol.jl")
include("ConvertToBig.jl")
include("GreensInt.jl")
mm=41 # number of values
Al=logspace(-3.5,-2,mm)
C_Guess=logspace(-20,-10,mm)
cl=Int(5)

acc=zeros(mm,mm)
#acc2=zeros(mm,1)
stab=zeros(mm,mm)



# Note: CL must be consistent.
# Note: choice of acc/stab means can be varied from poly to poly
# Note: value of c_1,c_2 can be varied from poly AND Alpha

Z=[];
for i=1:mm
    for j=1:mm
        #RUN2(NG,Al,n,d,cl,opt,plotopt,c_1,c_2)
        (acc[i,j],acc2,stab[i,j])=RUN2(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(5),Int(8),0,1/(C_Guess[j]),1)
    end
    println(i)
end




