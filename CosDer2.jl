function CosDer2(D,A,n,F)
    # first derivative of the cos expansion
    hc=BigFloat(1)/BigFloat((n-1));
    Coarse=collect(BigFloat,-D+hc:hc:D)
    Coarse=Coarse.-(5/90)
    LC=length(Coarse)
    hf=BigFloat(1)/(BigFloat((n-1))*BigFloat(F));
    Fine=collect(BigFloat,-D+hf:hf:D)
    Fine=Fine.-(5/900);
    F1=collect(BigFloat,A:hf:D-A);
    SB=length(F1);
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    #Modes=5;
    M=zeros(BigFloat,LF,Modes+1)
    B=zeros(BigFloat,SB,Modes+1);
        for j=1:Modes+1
            M[:,j]=-(j-1)^2*pi^2/D^2*cos.((j-1)*pi/D*Fine);
            B[:,j]=-(j-1)^2*pi^2/D^2*cos.((j-1)*pi/D*F1);
            
        end
    return M,B
end
