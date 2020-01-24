function expr=fmaker(d,n)
syms x
h=2/(n-1);
z=[-1:h:1]';
A=fliplr(vander(z));
[Q,R]=qr(A);
if d==0
    expr=1/R(1,1)*x^0;
else
    expr=x^d;
    for i=0:d-1
        qi=fmaker(i,n);
        expr=expr-qi*R(i+1,d+1);
    end
    expr=expr/R(d+1,d+1);
end
