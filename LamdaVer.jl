include("RUNStab.jl")
include("Spaces.jl")

Lam=logspace(-15,0,151);

AC=zeros(BigFloat,151,9);
Stab=zeros(BigFloat,151,9);

for j=1:9

for i=1:151
(Stab[i,j],AC[i,j],F1,F2,F3,F4)=RUNStab(100,.01,11,Int(j),5,27,0,Lam[i]);
println(i)
end
end


