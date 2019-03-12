

function linspace(a,b,n)
    z=zeros(BigFloat,n,1);
    h=(BigFloat(b-a)/BigFloat(n-1));
    for i=1:n
        z[i]=a+h*(i-1)
    end
return z
end

function logspace(a,b,n)
z=linspace(a,b,n)
q=10 .^ z
return q
end


