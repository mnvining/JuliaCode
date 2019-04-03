include("CosCont.jl")
include("StabBisect.jl")
using LinearAlgebra
using SymPy
using GenericSVD
using DelimitedFiles

K2=zeros(BigFloat,10,3)
Al=.01
for i = 1:10
    L2=StabBisect(i-1,Al)
    (a,s)=CosCont(i-1,Al,10^L2)
    K2[i,1]=L2
    K2[i,2]=s
    K2[i,3]=a
    println(i)
end

open("Al01.txt","w") do f
    writedlm(f,K2)
end



K3=zeros(BigFloat,10,3)
Al=.1
for i = 1:10
    L3=StabBisect(i-1,Al)
    (a,s)=CosCont(i-1,Al,10^L3)
    K3[i,1]=L3
    K3[i,2]=s
    K3[i,3]=a
    println(i)
end

open("Al1.txt","w") do f
    writedlm(f,K3)
end



