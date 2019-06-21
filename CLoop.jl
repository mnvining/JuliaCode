include("NCalc.jl")
CAM9=zeros(40,1)
SAM9=zeros(40,1)
UCM9=zeros(40,1)
USM9=zeros(40,1)
LL=collect(BigFloat,-10:29)


    for j=-10:29
  (CAM9[j+11],SAM9[j+11],UCM9[j+11],USM9[j+11])=NCalc(9,10 .^LL[j+11],FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm);
    println(j)
    end
