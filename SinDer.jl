function SinDer(D,A,n,F)
    # second derivative of the cos expansion
    hc=1/(n-1);
    Coarse=collect(-D+hc:hc:D)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(-D+hf:hf:D)
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,Modes,Modes)
        for j=1:Modes
            M[j,j]=(j)^2*pi^2/D^2;
        end
    return M
end
