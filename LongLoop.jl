include("CosContE.jl")
include("CosContLM.jl")

lam=logspace(-15,0,500);
mu=logspace(-15,0,500);
gam=logspace(-15,0,500);
global M= zeros(BigFloat,1,3);
	for i = 1:500
	    for j=1:500
	    	for k=1:500
		    (a,b,c)=CosContE(2,.001,lam[i],mu[j],gam[k])
                    if b<=1.0
                        nn=hcat(a,b,maximum(abs.(c)));
                        global M
                        M=hcat(M,nn);
                    end
                end
            end
        end



		   
