function CosMtx(D,A,n,F)
    # D = endpoint, period 2D
    # A = starting point
    # n = number of points, coarse grid
    # F = mult. factor coarse to fine grid


    hc=BigFloat(4)/BigFloat((35));
    Coarse=collect(BigFloat,-D+hc:hc:D)
    LC=length(Coarse)
    hf=BigFloat(hc)/(BigFloat(F));
    Fine=collect(BigFloat,-D+hc:hf:D+1e-10)
    F1=collect(BigFloat,A:hf:D-A+1e-10);
    SB=length(F1);
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    #Modes=5;
    M=zeros(BigFloat,LF,Modes+1)
    B=zeros(BigFloat,SB,Modes+1);
    C=zeros(BigFloat,LF,LF);

        for j=1:Modes+1
            M[:,j]=cos.((j-1)*Fine/D*pi);
            B[:,j]=cos.((j-1)*F1/D*pi);
        end
    for i=1:LF
        C[:,i]=cos.((i-1)*Fine/D*pi)
    end


    return M,B,C
end
