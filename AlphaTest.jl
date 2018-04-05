using PyPlot
using SymPy
using GenericSVD
include("TryinBig.jl")
include("All.jl")
include("dftmtx.jl")
mm=71 # number of points
Al=logspace(-5,2,mm)
#Al=0.01
#mm=5

    acc=zeros(mm,5,9)
    acc2=zeros(mm,5,9)
stab=zeros(mm,5,9)
for j=0:8
    
    for i=1:5
        println(i)
        for d=1:mm
            (acc[d,i,j+1],acc2[d,i,j+1],stab[d,i,j+1])=RUN2(Int(32),BigFloat(Al[d]),Int(9),Int(j),Int(i))         
        end
    end
end
        

#figure(1)
#loglog(Al,1-stab[:,1],Al,1-stab[:,2],Al,1-stab[:,3],Al,1-stab[:,4],Al,1-stab[:,5],Al,1-stab[:,6],Al,1-stab[:,7],Al,1-stab[:,8],Al,1-stab[:,9])

#figure(2)
#loglog(Al,acc[:,1],Al,acc[:,2],Al,acc[:,3],Al,acc[:,4],Al,acc[:,5],Al,acc[:,6],Al,acc[:,7],Al,acc[:,8],Al,acc[:,9])

#figure(3)
#loglog(Al,acc2[:,1],Al,acc2[:,2],Al,acc2[:,3],Al,acc2[:,4],Al,acc2[:,5],Al,acc2[:,6],Al,acc2[:,7],Al,acc2[:,8],Al,acc2[:,9])

#figure(1)
#loglog(Al,1-stab[:,1],Al,1-stab[:,2],Al,1-stab[:,3],Al,1-stab[:,4],Al,1-stab[:,5])

#figure(2)
#loglog(Al,acc[:,1],Al,acc[:,2],Al,acc[:,3],Al,acc[:,4],Al,acc[:,5])

#figure(3)
#loglog(Al,acc2[:,1],Al,acc2[:,2],Al,acc2[:,3],Al,acc2[:,4],Al,acc2[:,5])
