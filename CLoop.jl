include("NCalc.jl")
CAM3=zeros(40,1)
SAM3=zeros(40,1)
UCM3=zeros(40,1)
USM3=zeros(40,1)
LL=logspace(-10,29,40)


    for j=-10:29
  (CAM3[j+11],SAM3[j+11],UCM3[j+11],USM3[j+11])=NCalc(d,LL[j+11],FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm);
    println(j)
    end
