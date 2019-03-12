
using SymPy
using GenericSVD
using LinearAlgebra
function RUNStab(F,Al,n,d,cl,NoMo,opt)
    # uses NG grid points (int)
    # uses Al for alpha - should be a big float
    # n for number of coarse grid pts on -1,1 (int)
    # d for degree of gram poly (int)
    # cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    # opt: 0 for even, 1 for odd
    # lam - weight of h1,h2
    # mu - weight of u
    setprecision(200)
    


    include("ConvertToBig.jl")
    include("dftmtx.jl")
    include("Greens.jl")
    include("hsol.jl")
    include("coeffcalc.jl")
    include("Spaces.jl")
    include("All.jl")
    include("GreensInt.jl")
    include("FuncOut.jl")
    include("myhouse.jl")
    x=Sym("x")

    h1=2/(n-1);
    xcoarse1=collect(-1:h1:1-h1);
    #xcoarse2=collect(1+cl:h1:3+cl-h1);

    h2=h1/F;
    

    #h=xx[2]-xx[1]
    #x2=collect(6:h:1+2*cl-h1)
    #xall=collect(-1:h:1+2*cl-h1)
    xx=collect(-1:h2:1-h1);
    NG=length(xx);
    
    f(x)=fmaker(d,n)
    ef=ConvertToBig(evalasarray(f(x),xx))
    ef2=ConvertToBig(evalasarray((-1)^opt*f(x),xx))
    ef=vcat(ef,ef2)
    efff=ConvertToBig(evalasarray(f(x),xcoarse1))
    M=round(Int,round(Int,floor(NoMo/2))) # x2 because double modes
  
    ZZ=round(Int,((cl*NG/2+NG))) # x2 because double domain

    LD=round(Int,NG/(n-1))
    Modes=vcat(collect(1:M+1),collect(ZZ-M+1:ZZ))
    NMU=length(Modes);
    

    A=dftmtx(ZZ)
    A=A';
    A=A/sqrt(ZZ)

    A=A[:,Modes]
    
    ACont1=A[1:NG,:]
    #ACont2=A[round(Int,ZZ/2)+1:round(Int,ZZ/2)+NG,:]


    
    LL=round(Int,ZZ/2)
    fd=((cl+2)/2);
    K=zeros(BigFloat,ZZ)
    K=-pi*1im*vcat(collect(0:LL-1),0,collect(-(LL-1):-1))

    
    # derivative coefficients pertaining to the FC
    qq=vcat(K[1:M+1],K[(ZZ)-M+1:(ZZ)])/fd;
    D = Diagonal(I,(length(qq)))-Al*Diagonal{BigFloat}(qq.*qq)
    D_Dag=pinv(D,1e-40)

    GG(x)=totalGreen(d,n,Al,xcoarse1[end])
    
    (h1,h2)=hsol(xx,Al,xcoarse1[end])


    EGGf=ConvertToBig(evalasarray(GG(x),xx))
    Eh1f=ConvertToBig(evalasarray(h1(x),xx))
    Eh2f=ConvertToBig(evalasarray(h2(x),xx))
    U_G=ConvertToBig(evalasarray(GG(x),xcoarse1))


                     Ent_1=hcat(ACont1,zeros(NG,2))
                     #Ent_2=hcat(ACont2,zeros(NG,2))
                     Ent_3=hcat(ACont1*D_Dag,Eh1f,Eh2f)
                     #Ent_4=hcat(ACont2*D_Dag,(-1)^opt*Eh1f,(-1)^opt*Eh2f)


    #AugMat=vcat(Ent_1,Ent_2,Ent_3,Ent_4)
    AugMat=vcat(Ent_1,Ent_3);
                     #P=ACont1*D_Dag;
                     #P2=ACont1*D_Dag;
                     #P2[1:LD:end,:]=zeros(length(1:LD:NG),length(Modes));
    #US=P-P2;# this is the coarse grid eval of u I THINK

    #P=ACont2*D_Dag;
    #P2=ACont2*D_Dag;
    #P2[1:LD:end,:]=zeros(length(1:LD:NG),length(Modes));
    #US2=P-P2;


    #Ent_6=hcat(mu*US,zeros(NG,2));
    #Ent_7=hcat(mu*US2,zeros(NG,2));

    

                     


                     
                     

    #if lam>0
        

        #Ent_5=hcat(zeros(size(ACont1)),lam*Eh1f,lam*Eh2f)
        #AugMat=vcat(AugMat,Ent_5)
                     #end

    #AugMat=vcat(AugMat,Ent_6,Ent_7);


    if opt==0
        if mod(NMU,2)==0
            AugMat=AugMat[:,setdiff(1:end,vcat(2:2:M,M+1:2:end-2))]
        else
            AugMat=AugMat[:,setdiff(1:end,vcat(2:2:M+2,M+3:2:end-2))]
        end
    else
        if mod(NMU,2)==0
            AugMat=AugMat[:,setdiff(1:end,vcat(1:2:M,M+2:2:end-2))]
        else
            AugMat=AugMat[:,setdiff(1:end,vcat(1:2:M,M+3:2:end-2))]
        end
    end

 

    
#if lam >0
 #   AugVec=vcat(ef,EGGf,(-1)^opt*EGGf,0*Eh1f)
#else
    AugVec=vcat(ef[1:NG],EGGf)
    println(size(AugVec))
    println(size(AugMat))

                     #AugVec=vcat(AugVec,0*Eh1f,0*Eh1f)


    (U,S,V)=GenericSVD.svd(AugMat)
    S=Diagonal{BigFloat}(S)
    S_Dag=pinv(S,1e-40)
    
    Res=V*(S_Dag*(U'*AugVec))

    
    
    

    Cin=Res[1:end-2];

if opt == 0
    if mod(NMU,2)==0
       A=A[:,setdiff(1:end,vcat(2:2:M,M+1:2:end))]
    else
       A=A[:,setdiff(1:end,vcat(2:2:M+2,M+3:2:end))]
    end
else
    if mod(NMU,2)==0
        A=A[:,setdiff(1:end,vcat(1:2:M,M+2:2:end))]
    else
        A=A[:,setdiff(1:end,vcat(1:2:M,M+3:2:end))]
    end
end

if opt == 0
    if mod(NMU,2)==0
        ACont1=ACont1[:,setdiff(1:end,vcat(2:2:M,M+1:2:end))]
    else
        ACont1=ACont1[:,setdiff(1:end,vcat(2:2:M+2,M+3:2:end))]
    end
else
    if mod(NMU,2)==0
        ACont1=ACont1[:,setdiff(1:end,vcat(1:2:M,M+2:2:end))]
    else
        ACont1=ACont1[:,setdiff(1:end,vcat(1:2:M,M+3:2:end))]
    end
end
if opt == 0
    if mod(NMU,2)==0
        D_Dag=D_Dag[setdiff(1:end,vcat(2:2:M,M+1:2:end)),:]
    else
        D_Dag=D_Dag[setdiff(1:end,vcat(2:2:M+2,M+3:2:end)),:]
    end
else
    if mod(NMU,2)==0
        D_Dag=D_Dag[setdiff(1:end,vcat(1:2:M,M+2:2:end)),:]
    else
        D_Dag=D_Dag[setdiff(1:end,vcat(1:2:M,M+3:2:end)),:]
    end
end

if opt == 0
    if mod(NMU,2)==0
        D_Dag=D_Dag[:,setdiff(1:end,vcat(2:2:M,M+1:2:end))]
    else
        D_Dag=D_Dag[:,setdiff(1:end,vcat(2:2:M+2,M+3:2:end))]
    end
else
    if mod(NMU,2)==0
        D_Dag=D_Dag[:,setdiff(1:end,vcat(1:2:M,M+2:2:end))]
    else
        D_Dag=D_Dag[:,setdiff(1:end,vcat(1:2:M,M+3:2:end))]
    end
end


println(size(A))

    
    F1=real(A*Cin);
    FF1=real(ACont1*D_Dag*Cin)#+Res[end-1]*Eh1f+Res[end]*Eh2f);

acc=norm(F1[1:NG]-ef[1:NG])

#FF2=real(ACont1*D_Dag*Cin+Res[end-1]*Eh1f+Res[end]*Eh2f)



    
# coarse grid things
solcoarse=FF1[1:LD:end];

#ExAc=norm(FF2[1:100]-EGGf[1:100])

    
    stab=norm(solcoarse)/norm(efff)
    

    
    #NB=round(Int,floor(length(Eh1f)/2))

    return acc,stab,F1,EGGf,ef,efff#FF1,ef,Res[end-1:end],FF2,ExAc
end




