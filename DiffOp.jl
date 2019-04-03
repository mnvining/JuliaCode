function DiffOp(D,A,n,F,Al,opt)
    # opt 0 for sine
    # opt 1 for cosine
    hc=BigFloat(1)/BigFloat((n-1));
    Coarse=collect(BigFloat,-D+hc:hc:D)
    Coarse=Coarse.-(BigFloat(5)/BigFloat(90))
    LC=length(Coarse)
    hf=BigFloat(1)/(BigFloat((n-1))*BigFloat(F));
    Fine=collect(BigFloat,-D+hf:hf:D)
    Fine=Fine.-(BigFloat(5)/BigFloat(900));
    F1=collect(BigFloat,A:hf:D-A);
    SB=length(F1);
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);

    DO=zeros(BigFloat,Modes+opt,Modes+opt);
    if opt==1
    for i=1:Modes+1
        DO[i,i]=1+Al*(i-1)^2*pi^2/D^2
    end
    else
        for i =1:Modes
            DO[i,i]=1+Al*i^2*pi^2/D^2
        end
    end
    
    return DO
end


