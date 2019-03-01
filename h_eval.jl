using SymPy
function h_eval(z,Al,p)
    include("hsol.jl")
    x=Sym("x")
    # z is grid of points on which the function should be evaluated   
    # Al is alpha 
    # p is the last point on the coarse grid

    (h1,h2)=hsol(z,Al,p)

    h_1_e=evalasarray(h1(x),z)
    h_2_e=evalasarray(h2(x),z)

    return h_1_e,h_2_e
end

function h_on_grid(D,A,n,F,Al)
    hc=1/(n-1);
    Coarse=collect(-D+hc:hc:D)
    LC=length(Coarse)
    hf=hc/F;
    Fine=collect(-D+hf:hf:D)
    LF=length(Fine)
    Modes=round(Int,floor(LC/2)-1);
    M=zeros(BigFloat,LF,Modes)
    p=1;
    z=-1+hf:hf:1;
    return (H1,H2)=h_eval(z,Al,p)
end

