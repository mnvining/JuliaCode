function CosDer(D,A,n,F)
    # first derivative of the cos expansion
    hc=1/(n-1);
    Coarse=collect(-D+hc:hc:D)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(-D+hf:hf:D)
    F1=collect(A:hf:D-A);
    SB=length(F1)
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,LF,Modes+1)
    B=zeros(BigFloat,SB,Modes+1)
        for j=1:Modes+1
            M[:,j]=-(j-1)*pi/D*sin.((j-1)*pi/D*Fine);
            B[:,j]=-(j-1)*pi/D*sin.((j-1)*pi/D*F1);
        end
    return M,B
end
