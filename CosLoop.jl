include("CosCont2.jl")
include("CosContLM.jl")

lam=logspace(-15,-5,200);
L=zeros(200,3)
M=zeros(200,3)


for i=1:200
    (a,b,c,d,e)=CosCont2(2,.001,lam[i]);
    (f,g,h)=CosContLM(2,.001,0,lam[i])
    L[i,1]=b;
    L[i,2]=maximum(abs.(c));
    L[i,3]=a;
    M[i,1]=g;
    M[i,2]=maximum(abs.(h))
    M[i,3]=f;
end
 
 
