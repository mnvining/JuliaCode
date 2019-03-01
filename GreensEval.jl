function GreensEval(n,d,NC,Al)
    include("GreensInt.jl")
    include("All.jl")
    # n is no of pts eval'd on   
    # d is degree of poly  
    # NC is no of pts on original coarse gride for poly comp
    # Al is alpha

    hf=2/(n-1);
    
    hc=2/(NC-1);

    Fine=-1:hf:1-hc;


    
    expr=totalGreen(d,NC,Al,1-hc);
    
    
    
end
