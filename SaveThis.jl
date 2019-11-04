AlList=AllModule.logspace(-15,1,171);
w=AllModule.logspace(-2,15,181);
d=6;
WT=zeros(171,1)
MaxVec=zeros(171,1)
ACV=zeros(171,1)

for j=1:171
    Al=AlList[j];
    (IDOC,IDOS)=AllModule.IDOCalc(Al);
    global PP=zeros(1,4);
    global MaxCount=Inf;
    for i=1:121
        (CS,SS,CA,SA,MC,MS)=AllModule.TotalCalc(d,w[i],1e-7,AllModule.FCMat,AllModule.FCMatS,AllModule.C_C,AllModule.C_S,IDOC,IDOS,AllModule.m,AllModule.mm);
        if 1-CS>1e-18
            if MC<1e7
                IL=hcat(CS,CA,MC,i);
                global PP=vcat(PP,IL);
            end
        end
    end
    println(j)
    ACV[j]=PP[end,2];
    MaxVec[j]=PP[end,3];
    WT[j]=PP[end,end];
end
