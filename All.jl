function fmaker(d,n)
    x=Sym("x")
    h=BigFloat(2)/BigFloat(n-1);
    z=collect(BigFloat,-1:h:1);
    A=fliplr(myvander(z))
    (Q,R)=HRQR(A)
    #temp fix
    R[end,end]=R[end,end]/3
    if d==0
        expr=1/R[1,1]*x^0
        return expr
    else
        expr=x^d
        for i=0:d-1
            qi=fmaker(i,n)
            expr=expr-qi*R[i+1,d+1]
        end
        if d==9
            expr=expr*3
        end
        return expr/(R[d+1,d+1])
    end
end
