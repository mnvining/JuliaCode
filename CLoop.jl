include("NCalc.jl")
CAM=zeros(20,10)
SAM=zeros(20,10)
UCM=zeros(20,10)
USM=zeros(20,10)

for i=0:9
    for j=-1:18
  (CAM[j+2,i+1],SAM[j+2,i+1],UCM[j+2,i+1],USM[j+2,i+1])=NCalc(i,10.0^j,FCMat,FCMatS,C_C,C_S,IDOC,IDOS,m,mm);
    println(j)
    end
end
