include("NCalc.jl")
CAM5=zeros(40,1)
SAM5=zeros(40,1)
UCM5=zeros(40,1)
USM5=zeros(40,1)
LL=logspace(-10,29,40)


    for j=-10:29
  (CAM5[j+11],SAM5[j+11],UCM5[j+11],USM5[j+11])=NCalc(d,LL[j+11],FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm);
    println(j)
    end
