include("FuncOut.jl")
include("TryinBig.jl")
mm=121
Al=logspace(-4,2,mm)
z=zeros(1,2)
for i=1:mm
    (acc,acc2,stab,d)=RUN2(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(5),Int(8),0,1e5,1)
    if stab <=1
    L=hcat(Al[i],d)
        z=vcat(z,L)
    end
    
end
