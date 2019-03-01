function CosMtx(D,A,n,F)
    # D = endpoint, period 2D
    # A = starting point
    # n = number of points, coarse grid
    # F = mult. factor coarse to fine grid

    
    hc=1/(n-1);
    Coarse=collect(-D+hc:hc:D)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(-D+hf:hf:D)
    F1=collect(A:hf:D-A);
    SB=length(F1);
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,LF,Modes+1)
    B=zeros(BigFloat,SB,Modes+1);

        for j=1:Modes+1
            M[:,j]=cos.((j-1)*Fine/D*pi);
            B[:,j]=cos.((j-1)*F1/D*pi);
        end

    return M,B
end

    
    
