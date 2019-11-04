function IDOCalc(Al)
setprecision(230);
D=BigFloat(350)/BigFloat(100);
A=BigFloat(125)/BigFloat(100);
n=10;
F=10;
  (DOC)=DiffOp(D,A,n,F,Al,1)
  (DOS)=DiffOp(D,A,n,F,Al,0)
  IDOC=zeros(BigFloat,size(DOC))
  IDOS=zeros(BigFloat,size(DOS))
  for i=1:length(DOC[:,1])
    IDOC[i,i]=BigFloat(1)/DOC[i,i]
  end
  for i=1:length(DOS[:,1])
    IDOS[i,i]=BigFloat(1)/DOS[i,i]
  end
return(IDOC,IDOS)
end
