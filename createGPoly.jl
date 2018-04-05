function createGPoly(x,z,d)
    # d is degree
    # n is number of points


if d==0
    return 1/sqrt(sum(ones(size(x))));
else
    y(z)=z.^d;
    for i = 1:d-1
        ff(z)=createGPoly(x,z,i);
        y(z) = y(z)-sum(y(x).*ff(x))*ff(z);
    end
    return y(z)/sqrt(sum((f(x)).^2));
 
end




end
