Z=zeros(length(Al),1)
for i=1:length(Al)
    global k=0
    for j=1:length(w)
        if (1 .-cosstab[i,j])>(Al[i])
            global k=vcat(k,j)
        end
    end
    figure(1)
    loglog(w[k[2:end]],cosacc[i,k[2:end]])
    Z[i]=length(k[2:end])
    figure(2)
    loglog(w[k[2:end]],(1 .-cosstab[i,k[2:end]]))

end
