Al=BigFloat(0.01);
D=BigFloat(350)/BigFloat(100);
A=BigFloat(125)/BigFloat(100);
n=10;
F=10;

setprecision(230)
FCMat=zeros(BigFloat,31,10)
AMat=zeros(BigFloat,10,1)
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
