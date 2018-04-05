using PyPlot
using SymPy
using GenericSVD
include("TryinBig.jl")
include("All.jl")
include("dftmtx.jl")
mm=61 # number of points
Al=logspace(-4,2,mm)
#Al=0.01
#mm=5
#mm=1
cc=5
    acc=zeros(mm,8,5)
    acc2=zeros(mm,8,5)
stab=zeros(mm,8,5)
for c=1:5
    for j=0:7
        #for j=1:5
        for d=1:mm
            #for d=1:5
            (acc[d,j+1,c],acc2[d,j+1,c],stab[d,j+1,c])=RUN2(Int(32),BigFloat(Al[d]),Int(9),Int(j),Int(c),2,0)         
        end
    end
end


        

figure(1)
loglog(Al,1-stab[:,1,2],Al,1-stab[:,2,2],Al,1-stab[:,3,2],Al,1-stab[:,4,2],Al,1-stab[:,5,2],Al,1-stab[:,6,2],Al,1-stab[:,7,2],Al,1-stab[:,8,2])

figure(2)
loglog(Al,acc[:,1,2],Al,acc[:,2,2],Al,acc[:,3,2],Al,acc[:,4,2],Al,acc[:,5,2],Al,acc[:,6,2],Al,acc[:,7,2],Al,acc[:,8,2])

#figure(3)
#loglog(Al,acc2[:,1],Al,acc2[:,2],Al,acc2[:,3],Al,acc2[:,4],Al,acc2[:,5],Al,acc2[:,6],Al,acc2[:,7],Al,acc2[:,8],Al,acc2[:,9])

#figure(1)
#loglog(Al,1-stab[:,1],Al,1-stab[:,2],Al,1-stab[:,3],Al,1-stab[:,4],Al,1-stab[:,5])

#figure(2)
#loglog(Al,acc[:,1],Al,acc[:,2],Al,acc[:,3],Al,acc[:,4],Al,acc[:,5])

#figure(3)
#loglog(Al,acc2[:,1],Al,acc2[:,2],Al,acc2[:,3],Al,acc2[:,4],Al,acc2[:,5])
#figure()

#semilogy(0:7,1-stab[:,1],0:7,1-stab[:,2],0:7,1-stab[:,3],0:7,1-stab[:,4],0:7,1-stab[:,5])


#loglog(Al,1-stab[1,:],Al,1-stab[2,:],Al,1-stab[3,:],Al,1-stab[4,:],Al,1-stab[5,:],Al,1-stab[6,:],Al,1-stab[7,:],Al,1-stab[8,:])
