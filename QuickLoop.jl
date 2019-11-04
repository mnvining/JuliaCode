StoreMe=zeros(630,171)
for i=1:171
  Al=AlList[i]
  (IDOC,IDOS)=AllModule.IDOCalc(Al);
  (yc,ys,Gc,Gs)=AllModule.MP(2,w[WTI[i]],AllModule.FCMat,AllModule.FCMatS,AllModule.C_C,AllModule.C_S,IDOC,IDOS,1e-20);
  StoreMe[:,i]=AllModule.m*yc;
  println(i)
end
