function SinMtx(D,A,n,F)
    # D = endpoint, period 2D
    # A = starting point
    # n = number of points, coarse grid
    # F = mult. factor coarse to fine grid


    hc=BigFloat(4)/BigFloat((35));
    Coarse=collect(BigFloat,-D+hc:hc:D)
    #Coarse=Coarse.-(5/90)
    LC=length(Coarse)
    hf=BigFloat(hc)/(BigFloat(F));
    Fine=collect(BigFloat,-D+hc:hf:D+1e-10)
    #Fine=Fine.-5/900;
    LF=length(Fine)
    F1=collect(BigFloat,A:hf:D-A+1e-10);
    SB=length(F1)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,LF,Modes)
    B=zeros(BigFloat,SB,Modes)
    C=zeros(BigFloat,LF,LF)


        for j=1:Modes
            M[:,j]=sin.((j)*Fine/D*pi);
            B[:,j]=sin.((j)*F1/D*pi);
        end
    for i=1:LF
        C[:,i]=sin.(i*Fine/D*pi)
        end

    return M,B,C
end
