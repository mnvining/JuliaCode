using PyPlot
using SymPy
using GenericSVD
include("TryinBig.jl")
include("All.jl")
include("dftmtx.jl")
include("coeffcalc.jl")
include("hsol.jl")
include("ConvertToBig.jl")
include("GreensInt.jl")
mm=41 # number of points
#Al=logspace(-3,2,mm)
C_Guess=0.00001:0.00001:1
cl=Int(5)
Al=3.1622776601683795

#acc=zeros(mm,1)
#acc2=zeros(mm,1)
#stab=zeros(mm,1)



# Note: CL must be consistent.
# Note: choice of acc/stab means can be varied from poly to poly
# Note: value of c_1,c_2 can be varied from poly AND Alpha
j=1;
Z=[];
for i=1:length(C_Guess)
    #for d=1:mm
        #RUN2(NG,Al,n,d,cl,opt,plotopt,c_1,c_2)
    (acc,acc2,stab)=RUN2(Int(32),BigFloat(Al),Int(9),Int(6),Int(5),Int(8),0,1/C_Guess[i],C_Guess[i])
    if stab<=1
        Z=vcat(Z,hcat(C_Guess[i],stab))
    end

    

    #end

    println(i)
end

Res=hcat(Al,acc,stab)



        

#figure(1)
#loglog(Al,1-stab[:,1],Al,1-stab[:,2],Al,1-stab[:,3],Al,1-stab[:,4],Al,1-stab[:,5],Al,1-stab[:,6],Al,1-stab[:,7],Al,1-stab[:,8])

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
#figure()

#semilogy(0:7,1-stab[:,1],0:7,1-stab[:,2],0:7,1-stab[:,3],0:7,1-stab[:,4],0:7,1-stab[:,5])


#loglog(Al,1-stab[1,:],Al,1-stab[2,:],Al,1-stab[3,:],Al,1-stab[4,:],Al,1-stab[5,:],Al,1-stab[6,:],Al,1-stab[7,:],Al,1-stab[8,:])
