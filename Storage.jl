setprecision(230)
Al=BigFloat(1)/BigFloat(1000);
# al = 1e-3
D=BigFloat(4)
A=BigFloat(52)/BigFloat(35);
n=10;
F=10;


FCMat=zeros(BigFloat,35,10)
AMat=zeros(BigFloat,10,1)
FCMatS=zeros(BigFloat,34,10)
AMatS=zeros(BigFloat,10,1)
for i=0:9
  (a,b,c)=CosCont(i)
  FCMat[:,i+1]=b;
  AMat[i+1]=a;
end


DOC=DiffOp(D,A,n,F,Al,1);
IDOC=zeros(BigFloat,size(DOC))
for i=1:length(DOC[1,:])
    IDOC[i,i]=1/DOC[i,i];
end

for i=0:9
  (a,b,c)=SinCont(i)
  FCMatS[:,i+1]=b;
  AMatS[i+1]=a;
end


DOS=DiffOp(D,A,n,F,Al,0);
IDOS=zeros(BigFloat,size(DOS))
for i=1:length(DOS[1,:])
    IDOS[i,i]=1/DOS[i,i];
end

(m,C_C,k)=CosMtx(D,A,n,F);
(mm,C_S,kk)=SinMtx(D,A,n,F);
