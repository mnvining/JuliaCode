include("CosCont.jl")
include("SinCont.jl")
include("DiffOp.jl")
include("CosMtx.jl")
include("SinMtx.jl")
using DelimitedFiles

Al=BigFloat(0.01);
D=BigFloat(350)/BigFloat(100);
A=BigFloat(125)/BigFloat(100);
n=10;
F=10;

setprecision(230)
FCMat=zeros(BigFloat,31,10)
AMat=zeros(BigFloat,10,1)
FCMatS=zeros(BigFloat,30,10)
AMatS=zeros(BigFloat,10,1)
for i=0:9
  (a,b,c)=CosCont(i)
  FCMat[:,i+1]=b;
  AMat[i+1]=a;
end

open("FcCos.txt","w") do k
    writedlm(k,FCMat)
end
open("CosAcc.txt","w") do s
    writedlm(s,AMat)
end

DOC=DiffOp(D,A,n,F,Al,1);
IDOC=zeros(BigFloat,size(DOC))
for i=1:length(DOC[1,:])
    IDOC[i,i]=1/DOC[i,i];
end
open("CosDiffOp.txt","w") do l
    writedlm(l,IDOC)
end
for i=0:9
  (a,b,c)=SinCont(i)
  FCMatS[:,i+1]=b;
  AMatS[i+1]=a;
end
open("FcSin.txt","w") do k
    writedlm(k,FCMatS)
end
open("SinAcc.txt","w") do s
    writedlm(s,AMatS)
end

DOS=DiffOp(D,A,n,F,Al,0);
IDOS=zeros(BigFloat,size(DOS))
for i=1:length(DOS[1,:])
    IDOS[i,i]=1/DOS[i,i];
end
open("SinDiffOp.txt","w") do l
    writedlm(l,IDOS)
end

(m,C_C,n)=CosMtx(D,A,n,F);
(mm,C_S,nn)=SinMtx(D,A,n,F);
