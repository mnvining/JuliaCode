include("FuncOut.jl")
include("TryinBig.jl")
mm=91
Al=logspace(-6,3,mm)
z=zeros(1,2)
for i=1:length(Al)
    (acc,acc2,stab,d)=RUN2(Int(32),BigFloat(Al[i]),Int(9),Int(6),Int(5),Int(8),0,1e10,1)
    if stab <=1
    L=hcat(Al[i],d)
        z=vcat(z,L)
    end

end
