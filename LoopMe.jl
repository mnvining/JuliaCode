include("RUNStab.jl")
include("Spaces.jl")
include("dftmtx.jl")
include("All.jl")
SMat=zeros(700,10);
LME=[1e-7 0;2e-8 0;3e-5 3e-5;0 1e-4;0 3.5e-4;0 1e-3;0 3e-3;0 7e-1;0 1.5;0 4];
for i=0:9
(s,a,F1,F2,F3,F4)=RUNStab(100,.001,11,i,5,27,1,LME[i+1,1],LME[i+1,2])
SMat[:,i+1]=F1;
end

