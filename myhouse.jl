using LinearAlgebra

function myhouse(x)
    m=length(x);
    s=x[2:m]'*x[2:m]
    v=vcat(1,x[2:m])
    if (s==0) & (x[1]>=0)
        b=0
    elseif (s==0) & (x[1]<0)
        b=-2 
    else  
        mu=sqrt(x[1]^2+s)
        if (x[1]<=0)
            v[1]=x[1]-mu   
        else   
            v[1]=-s/(x[1]+mu)    
        end 
        b=2*v[1]^2/(s+v[1]^2) 
        v=v/v[1]  
    end
    return v,b
end

function HRQR(A)
    (m,n)=size(A);
    A2=copy(A)
    Pr=zeros(BigFloat,size(A))
    Q1=zeros(BigFloat,(m,m))+Diagonal{BigFloat}(I,m);
    Q2=zeros(BigFloat,(m,m))+Diagonal{BigFloat}(I,m);
    for j=1:n
        (v,B)=myhouse(A2[j:m,j])
        s=length(v)
        II=zeros(BigFloat,(s,s))+Diagonal{BigFloat}(I,s);
        A2[j:m,j:n]=(II-B*(v*v'))*A2[j:m,j:n];
        if j < m
            A2[j+1:m,j]=v[2:m-j+1]
        end
        # Q matrix
        vj=vcat(zeros(j-1,1),v)
        Q1=Q1*(Diagonal{BigFloat}(I,m)-B*(vj*vj'));
    end
    for i = 1:m
        Pr[i,i:n]=A2[i,i:n]
    end
    
    Q2=Q2[:,1:n];
    for j=n:-1:1
        V=zeros(BigFloat,m,1);
        V[j:m]=vcat(1,A2[j+1:m,j]);
        B=2/(1+norm(A2[j+1:m,j],2)^2);
        Q2[j:m,j:n]=Q2[j:m,j:n]-B*(V[j:m]*V[j:m]')*Q2[j:m,j:n]
    end                             
    R=Pr;
    
    return Q1[:,1:n],R[1:n,1:n]

end


function BackSolve(R,y)
    # square R
    (m,n)=size(R);
    x=zeros(BigFloat,m,1);
    for i=m:-1:1
        x[i]=copy(y[i])
        for j=i+1:m
            x[i]=x[i]-R[i,j]*x[j]
        end
        x[i]=x[i]/R[i,i]
    end
    return x
end



function SolveViaQR(A,b)
    (Q,R)=HRQR(A)
    yy=Q'*b
    x=BackSolve(R,yy)
    return x,yy
end



        

