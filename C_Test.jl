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
mm=61 # number of values
Al=logspace(-4,2,mm)
C_Guess=logspace(1,7,mm)
cl=Int(5)

acc=zeros(mm,mm)
#acc2=zeros(mm,1)
stab=zeros(mm,mm)



# Note: CL must be consistent.
# Note: choice of acc/stab means can be varied from poly to poly
# Note: value of c_1,c_2 can be varied from poly AND Alpha

Z=zeros(1,3)
for i=1:mm
    for j=1:mm
        #RUN2(NG,Al,n,d,cl,opt,plotopt,c_1,c_2)
        (acc[i,j],acc2,stab[i,j],d)=RUN2(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(5),Int(8),0,(C_Guess[j],1)
if stab[i,j]<=1
La=hcat(Al[i],C_Guess[j],stab[i,j])
Z=vcat(Z,La)
end
end
    
end




