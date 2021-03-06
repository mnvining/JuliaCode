using PyPlot
using SymPy
using GenericSVD
function RUN2(NG,Al,n,d,cl,opt,plotopt,c_1,c_2)
    # uses NG grid points (int)
    # uses Al for alpha - should be a big float
    # n for number of coarse grid pts on -1,1 (int)
    # d for degree of gram poly (int)
    # cl continuation length (extended domain is [-1,1+cl)- cl must be type Int
    # opt: 1 for func, 2 for 2nd der, 3, for sobolev
    # c_1 - coeff of f c_2 coeff d^2 f
    
    


include("ConvertToBig.jl")
include("dftmtx.jl")
include("Greens.jl")
include("hsol.jl")
include("coeffcalc.jl")
include("All.jl")
    include("GreensInt.jl")
    include("FuncOut.jl")
x=Sym("x")

xx=linspace(-1,1,NG+1)
    xx=xx[1:end-1]
    h=xx[2]-xx[1]
    x2=collect(-1:h:1+cl-h)
    xcoarse1=collect(linspace(-1,1,n))
    xcoarse=xcoarse1[1:end-1]
    hh=xcoarse[2]-xcoarse[1]
    xcoarseext=collect(-1:hh:1+cl-hh)
f(x)=fmaker(d,n)
    ef=ConvertToBig(evalasarray(f(x),xx))
    efff=ConvertToBig(evalasarray(f(x),xcoarse))
M=round(Int,floor(length(xcoarseext)/2))
  
ZZ=Int(cl*NG/2+NG);
    Modes=vcat(collect(1:M+1),collect(ZZ-M+1:ZZ))
    #println(length(Modes))

A=createFWDcont(NG,Modes,cl)
B=createBackDoubleDomain(NG,Modes,cl)

    GG(x)=totalGreen(d,n,Al)
    d2GG(x)=diff(GG(x),x,2)
    GPrime(x)=diff(GG(x),x,1)
    
    (h1,h2)=hsol(xx,Al)
    h1Prime(x)=diff(h1(x),x,1)
    h2Prime(x)=diff(h2(x),x,1)
    d2h1(x)=diff(h1(x),x,2)
    d2h2(x)=diff(h2(x),x,2)

    
    EGG=ConvertToBig(evalasarray(GG(x),xcoarse1))
    Eh1=ConvertToBig(evalasarray(h1(x),xcoarse1))
    Eh2=ConvertToBig(evalasarray(h2(x),xcoarse1))
    eh1pr=ConvertToBig(evalasarray(h1Prime(x),xcoarse1))
    eh2pr=ConvertToBig(evalasarray(h2Prime(x),xcoarse1))
    Gpr=ConvertToBig(evalasarray(GPrime(x),xcoarse1))
    EG2=ConvertToBig(evalasarray(d2GG(x),xcoarse1))
    EH1=ConvertToBig(evalasarray(d2h1(x),xcoarse1))
    EH2=ConvertToBig(evalasarray(d2h2(x),xcoarse1))
    
if Al<-10^(-3.5) || Al>10^(-1.5)
    if opt==1
        Z=[Eh1'*Eh1 Eh2'*Eh1; Eh1'*Eh2 Eh2'*Eh2]
zz=[-EGG'*Eh1; -EGG'*Eh2]
    elseif opt==2
        Z=[EH1'*EH1 EH2'*EH1; EH1'*EH2 EH2'*EH2]
zz=[-EG2'*EH1; -EG2'*EH2]
    elseif opt==3
        Z=[(Eh1'*Eh1+eh1pr'*eh1pr) (Eh1'*Eh2+eh1pr'*eh2pr); (Eh1'*Eh2+eh1pr'*eh2pr) (Eh2'*Eh2+eh2pr'*eh2pr)]
        zz=[-EGG'*Eh1-Gpr'*eh1pr; -EGG'*Eh2-Gpr'*eh2pr]
    elseif opt==4
        Z=[(Eh1'*Eh1+eh1pr'*eh1pr+EH1'*EH1) (Eh1'*Eh2+eh1pr'*eh2pr+EH1'*EH2); (Eh1'*Eh2+eh1pr'*eh2pr+EH1'*EH2) (Eh2'*Eh2+eh2pr'*eh2pr+EH2'*EH2)]
        zz=[-EGG'*Eh1-Gpr'*eh1pr-EG2'*EH1; -EGG'*Eh2-Gpr'*eh2pr-EG2'*EH1]
    else
        Z=[(c_1^2*Eh1'*Eh1+c_1*c_2*Eh1'*EH1+c_1*c_2*EH1'*Eh1+c_1*c_2*Eh1'*EH1+c_2^2*EH1'*EH1) (c_1^2*Eh2'*Eh1+c_1*c_2*Eh2'*EH1+c_1*c_2*EH2'*EH1+c_1*c_2*Eh2'*EH1+c_2^2*EH2'*EH1);(c_1^2*Eh1'*Eh2+c_1*c_2*Eh1'*EH2+c_1*c_2*EH1'*Eh2+c_1*c_2*Eh1'*EH2+c_2^2*EH1'*EH2) (c_1^2*Eh2'*Eh2+c_1*c_2*Eh2'*EH2+c_1*c_2*EH2'*EH2+c_1*c_2*Eh2'*EH2+c_2^2*EH2'*EH2)]
        zz=[-c_1^2*EGG'Eh1-c_1*c_2*EGG'*EH1-c_1*c_2*EG2'*Eh1-c_1*c_2*EGG'*Eh1-c_2^2*EG2'*EH1;-c_1^2*EGG'Eh2-c_1*c_2*EGG'*EH2-c_1*c_2*EG2'*Eh2-c_1*c_2*EGG'*Eh2-c_2^2*EG2'*EH2]
        end
    cobb=Z\zz

    

    EGGf=ConvertToBig(evalasarray(GG(x),xx))
    Eh1f=ConvertToBig(evalasarray(h1(x),xx))
    Eh2f=ConvertToBig(evalasarray(h2(x),xx))


    BB=EGGf+cobb[1]*Eh1f+cobb[2]*Eh2f
else
    AlVec=10.^[-3.4,-1.5,-1.45, -1.4, -1.35, -1.3, -1.25, -1.2, -1.15, -1.1, -1.05, -1.0, -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0];

    K=zeros(NG,length(AlVec));
    for jj=1:length(AlVec)
        K[:,jj]=FuncOut(NG,AlVec[jj],n,d,cl,opt,plotopt,c_1,c_2)
    end


    BB=1.137024986*K[:,1]-92.95971745*K[:,2]+117.6944579*K[:,3]-31.53613496*K[:,4]-0.1546041968e-24*K[:,50]-0.5222920362e-36*K[:,70]+Al^3*(0.7127261081e-1*K[:,1]-3.505334485*K[:,2]+4.354215980*K[:,3]-1.166708656*K[:,4]-0.5719726113e-26*K[:,50]-0.1932267984e-37*K[:,70]+.3126186421*K[:,5]-0.8376591268e-1*K[:,6]+0.2244500866e-1*K[:,7]-0.6014121944e-2*K[:,8]+0.1611479118e-2*K[:,9]-0.4317945283e-3*K[:,10]+0.1156989951e-3*K[:,11]-0.3100145231e-4*K[:,12]+0.8306814112e-5*K[:,13]-0.2225804133e-5*K[:,14]+0.5964024199e-6*K[:,15]-0.1598055468e-6*K[:,16]+0.4281976721e-7*K[:,17]-0.1147352204e-7*K[:,18]+0.3074320966e-8*K[:,19]-0.8237618201e-9*K[:,20]+0.2207263144e-9*K[:,21]-0.5914343770e-10*K[:,22]+0.1584743637e-10*K[:,23]-0.4246307778e-11*K[:,24]+0.1137794740e-11*K[:,25]-0.3048711817e-12*K[:,26]+0.8168998693e-13*K[:,27]-0.2188876603e-13*K[:,28]+0.5865077180e-14*K[:,29]-0.1571542694e-14*K[:,30]+0.4210935957e-15*K[:,31]-0.1128316889e-15*K[:,32]+0.3023315992e-16*K[:,33]-0.8100950786e-17*K[:,34]+0.2170643221e-17*K[:,35]-0.5816220982e-18*K[:,36]+0.1558451715e-18*K[:,37]-0.4175858785e-19*K[:,38]+0.1118917989e-19*K[:,39]-0.2998131716e-20*K[:,40]+0.8033469720e-21*K[:,41]-0.2152561724e-21*K[:,42]+0.5767771756e-22*K[:,43]-0.1545469784e-22*K[:,44]+0.4141073806e-23*K[:,45]-0.1109597382e-23*K[:,46]+0.2973157225e-24*K[:,47]-0.7966550773e-25*K[:,48]+0.2134630846e-25*K[:,49]+0.1532595993e-26*K[:,51]-0.4106578587e-27*K[:,52]+0.1100354416e-27*K[:,53]-0.2948390771e-28*K[:,54]+0.7900189262e-29*K[:,55]-0.2116849333e-29*K[:,56]+0.5672080692e-30*K[:,57]-0.1519829441e-30*K[:,58]+0.4072370713e-31*K[:,59]-0.1091188444e-31*K[:,60]+0.2923830623e-32*K[:,61]-0.7834380542e-33*K[:,62]+0.2099215938e-33*K[:,63]-0.5624832102e-34*K[:,64]+0.1507169028e-34*K[:,65]-0.4038440087e-35*K[:,66]+0.1082070071e-35*K[:,67]-0.2898401976e-36*K[:,68]+0.7729071937e-37*K[:,69]+0.3220446640e-38*K[:,71])+Al^2*(.7269806303*K[:,1]-35.75441175*K[:,2]+44.41300300*K[:,3]-11.90042829*K[:,4]-0.5834120636e-25*K[:,50]-0.1970913344e-36*K[:,70]+3.188710149*K[:,5]-.8544123093*K[:,6]+.2289390883*K[:,7]-0.6134404382e-1*K[:,8]+0.1643708700e-1*K[:,9]-0.4404304188e-2*K[:,10]+0.1180129751e-2*K[:,11]-0.3162148136e-3*K[:,12]+0.8472950394e-4*K[:,13]-0.2270320216e-4*K[:,14]+0.6083304683e-5*K[:,15]-0.1630016577e-5*K[:,16]+0.4367616255e-6*K[:,17]-0.1170299248e-6*K[:,18]+0.3135807385e-7*K[:,19]-0.8402370565e-8*K[:,20]+0.2251408407e-8*K[:,21]-0.6032630646e-9*K[:,22]+0.1616438510e-9*K[:,23]-0.4331233933e-10*K[:,24]+0.1160550635e-10*K[:,25]-0.3109686053e-11*K[:,26]+0.8332378667e-12*K[:,27]-0.2232654135e-12*K[:,28]+0.5982378724e-13*K[:,29]-0.1602973548e-13*K[:,30]+0.4295154676e-14*K[:,31]-0.1150883227e-14*K[:,32]+0.3083782312e-15*K[:,33]-0.8262969802e-16*K[:,34]+0.2214056086e-16*K[:,35]-0.5932545401e-17*K[:,36]+0.1589620749e-17*K[:,37]-0.4259375961e-18*K[:,38]+0.1141296349e-18*K[:,39]-0.3058094350e-19*K[:,40]+0.8194139115e-20*K[:,41]-0.2195612958e-20*K[:,42]+0.5883127191e-21*K[:,43]-0.1576379180e-21*K[:,44]+0.4223895282e-22*K[:,45]-0.1131789330e-22*K[:,46]+0.3032620369e-23*K[:,47]-0.8125881788e-24*K[:,48]+0.2177323463e-24*K[:,49]+0.1563247913e-25*K[:,51]-0.4188710158e-26*K[:,52]+0.1122361504e-26*K[:,53]-0.3007358587e-27*K[:,54]+0.8058193047e-28*K[:,55]-0.2159186319e-28*K[:,56]+0.5785522306e-29*K[:,57]-0.1550226030e-29*K[:,58]+0.4153818127e-30*K[:,59]-0.1113012213e-30*K[:,60]+0.2982307236e-31*K[:,61]-0.7991068153e-32*K[:,62]+0.2141200257e-32*K[:,63]-0.5737328744e-33*K[:,64]+0.1537312408e-33*K[:,65]-0.4119208889e-34*K[:,66]+0.1103711473e-34*K[:,67]-0.2956370016e-35*K[:,68]+0.7883653376e-36*K[:,69]+0.3284855573e-37*K[:,71])+Al*(1.688124229*K[:,1]-108.3844267*K[:,2]+135.2854905*K[:,3]-36.24963793*K[:,4]-0.1777118903e-24*K[:,50]-0.6003556627e-36*K[:,70]+9.713061209*K[:,5]-2.602606907*K[:,6]+.6973664189*K[:,7]-.1868587688*K[:,8]+0.5006865619e-1*K[:,9]-0.1341585599e-1*K[:,10]+0.3594767779e-2*K[:,11]-0.9632151234e-3*K[:,12]+0.2580927145e-3*K[:,13]-0.6915573441e-4*K[:,14]+0.1853022319e-4*K[:,15]-0.4965158339e-5*K[:,16]+0.1330410167e-5*K[:,17]-0.3564823299e-6*K[:,18]+0.9551915241e-7*K[:,19]-0.2559427975e-7*K[:,20]+0.6857966590e-8*K[:,21]-0.1837586609e-8*K[:,22]+0.4923798480e-9*K[:,23]-0.1319327826e-9*K[:,24]+0.3535128257e-10*K[:,25]-0.9472347615e-11*K[:,26]+0.2538107894e-11*K[:,27]-0.6800839605e-12*K[:,28]+0.1822279480e-12*K[:,29]-0.4882783150e-13*K[:,30]+0.1308337802e-13*K[:,31]-0.3505680574e-14*K[:,32]+0.9393442789e-15*K[:,33]-0.2516965409e-15*K[:,34]+0.6744188488e-16*K[:,35]-0.1807099859e-16*K[:,36]+0.4842109479e-17*K[:,37]-0.1297439324e-17*K[:,38]+0.3476478192e-18*K[:,39]-0.9315195241e-19*K[:,40]+0.2495999042e-19*K[:,41]-0.6688009277e-20*K[:,42]+0.1792046685e-20*K[:,43]-0.4801774619e-21*K[:,44]+0.1286631632e-21*K[:,45]-0.3447519066e-22*K[:,46]+0.9237599497e-23*K[:,47]-0.2475207325e-23*K[:,48]+0.6632298039e-24*K[:,49]+0.4761775750e-25*K[:,51]-0.1275913967e-25*K[:,52]+0.3418801170e-26*K[:,53]-0.9160650127e-27*K[:,54]+0.2454588804e-27*K[:,55]-0.6577050877e-28*K[:,56]+0.1762315471e-28*K[:,57]-0.4722110073e-29*K[:,58]+0.1265285581e-29*K[:,59]-0.3390322495e-30*K[:,60]+0.9084341746e-31*K[:,61]-0.2434142034e-31*K[:,62]+0.6522263920e-32*K[:,63]-0.1747635334e-32*K[:,64]+0.4682774169e-33*K[:,65]-0.1254743335e-33*K[:,66]+0.3361991711e-34*K[:,67]-0.9005334941e-35*K[:,68]+0.2401422651e-35*K[:,69]+0.1000592771e-36*K[:,71])+8.450081895*K[:,5]-2.264192620*K[:,6]+.6066885840*K[:,7]-.1625617161*K[:,8]+0.4355828056e-1*K[:,9]-0.1167140610e-1*K[:,10]+0.3127343839e-2*K[:,11]-0.8379692561e-3*K[:,12]+0.2245331854e-3*K[:,13]-0.6016348571e-4*K[:,14]+0.1612075741e-4*K[:,15]-0.4319543930e-5*K[:,16]+0.1157418308e-5*K[:,17]-0.3101293008e-6*K[:,18]+0.8309889571e-7*K[:,19]-0.2226628200e-7*K[:,20]+0.5966232279e-8*K[:,21]-0.1598647121e-8*K[:,22]+0.4283562051e-9*K[:,23]-0.1147776992e-9*K[:,24]+0.3075459182e-10*K[:,25]-0.8240668041e-11*K[:,26]+0.2208080347e-11*K[:,27]-0.5916533457e-12*K[:,28]+0.1585330362e-12*K[:,29]-0.4247879902e-13*K[:,30]+0.1138215989e-13*K[:,31]-0.3049840551e-14*K[:,32]+0.8172023128e-15*K[:,33]-0.2189686998e-15*K[:,34]+0.5867248627e-16*K[:,35]-0.1572124531e-16*K[:,36]+0.4212494986e-17*K[:,37]-0.1128734630e-17*K[:,38]+0.3024435325e-18*K[:,39]-0.8103950028e-19*K[:,40]+0.2171446865e-19*K[:,41]-0.5818374340e-20*K[:,42]+0.1559028706e-20*K[:,43]-0.4177404827e-21*K[:,44]+0.1119332250e-21*K[:,45]-0.2999241724e-22*K[:,46]+0.8036443978e-23*K[:,47]-0.2153358674e-23*K[:,48]+0.5769907177e-24*K[:,49]+0.4142606969e-25*K[:,51]-0.1110008192e-25*K[:,52]+0.2974257986e-26*K[:,53]-0.7969500255e-27*K[:,54]+0.2135421157e-27*K[:,55]-0.5721843746e-28*K[:,56]+0.1533163411e-28*K[:,57]-0.4108098978e-29*K[:,58]+0.1100761804e-29*K[:,59]-0.2949482364e-30*K[:,60]+0.7903114174e-31*K[:,61]-0.2117633061e-31*K[:,62]+0.5674180681e-32*K[:,63]-0.1520392117e-32*K[:,64]+0.4073877882e-33*K[:,65]-0.1091590356e-33*K[:,66]+0.2924835402e-34*K[:,67]-0.7834380542e-35*K[:,68]+0.2089168145e-35*K[:,69]+0.8704867269e-37*K[:,71]

    


end
    

D=createId2Derivative(NG,M,Al,cl)
D_Dag=pinv(D,1e-14)

AugMat=vcat(A,A*D_Dag')

AugVec=vcat(ef,BB)



Cin=AugMat\AugVec


    FF1ex=real(B*D_Dag*Cin)
    FF1=real(A*D_Dag*Cin)
    FF2=real(A*Cin)
    FF2ex=real(B*Cin)
acc=norm(FF1-BB)
    acc2=norm(FF2-ef)
    # coarse grid things
    LD=round(Int,NG/(n-1))
    solcoarse=FF1[1:LD:end];

    
    stab=norm(solcoarse)/norm(efff)
    #figure(d+1)
    if plotopt==1
        figure()

        plot(xx,FF1,xx,FF2,xx,ef,xx,BB,xx,EGGf)
        legend(("FullOp","Fc","func","Bf","Greens"))
    elseif plotopt==2
        plot(xx,ef,xx,BB,xx,EGGf)
        legend(("Func","Bf","Greens"))
    end
    
    

    
    NB=round(Int,floor(length(BB)/2))

    return acc,acc2,stab,BB
end




