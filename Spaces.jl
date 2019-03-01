

function linspace(a,b,n)
h=(b-a)/(n-1);
z=collect(a:h:b)
return z
end

function logspace(a,b,n)
z=linspace(a,b,n)
q=10 .^ z
return q
end


