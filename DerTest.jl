using GenericSVD
include("MP.jl")
include("NCalc.jl")
include("DiffOp.jl")
include("IDOCalc.jl")

#This is to check the inner product of f_C and f_C''


Al=1e-6;
setprecision(230);
D=BigFloat(350)/BigFloat(100);
A=BigFloat(125)/BigFloat(100);
n=10;
F=10;
D_O=DiffOp(D,A,n,F,Al,1);
(m,n)=size(D_O)
D_O=D_O-Diagonal{BigFloat}(I,m)

(IDOC,IDOS)=IDOCalc(Al)
DComp=zeros(BigFloat,size(IDOC))

(m,n)=size(IDOC)
for i=1:m
    DComp[i,i]=1-Al*(i-1)^2*pi^2/D^2
end
u=zeros(10,1)
q=zeros(10,1)
ll=zeros(10,1)
ww=zeros(10,1)
oo=zeros(10,1)
u1=zeros(10,1)
q1=zeros(10,1)
ll1=zeros(10,1)
ww1=zeros(10,1)
oo1=zeros(10,1)
for i = 1:10
    FC=C_C*FCMat[:,i];
    F1=C_C*D_O*FCMat[:,i];
    F2=C_C*DComp*FCMat[:,i];
    u[i]=norm(F1-F2);
    l1=FC[1:10:end];
    l2=F1[1:10:end];
    q[i]=l1'*l2;
    ll[i]=norm(F1[1:10:end])
    ww[i]=norm(FC[1:10:end])
    oo[i]=ww[i]+Al*ww[i]+Al^2*ll[i];
end
for i = 1:10
    (yc,ys,Gc,Gs)=MP(i-1,1e8,FCMat,FCMatS,C_C,C_S,IDOC,IDOS)
    FC=C_C*yc;
    F1=C_C*D_O*yc;
    F2=C_C*DComp*yc;
    u1[i]=norm(F1-F2);
    l1=FC[1:10:end];
    l2=F1[1:10:end];
    q1[i]=l1'*l2;
    ll1[i]=norm(F1[1:10:end])
    ww1[i]=norm(FC[1:10:end])
    oo1[i]=ww1[i]+Al*ww1[i]+Al^2*ll1[i];
end
