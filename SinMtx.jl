function SinMtx(D,A,n,F)
    # D = endpoint, period 2D
    # A = starting point
    # n = number of points, coarse grid
    # F = mult. factor coarse to fine grid

    
    hc=1/(n-1);
    Coarse=collect(-D+hc:hc:D)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(-D+hf:hf:D)
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,LF,Modes)

    for i=1:LF
        for j=1:Modes
            M[i,j]=sin((j)*Fine[i]/D*pi);
        end
    end
    return M
end
