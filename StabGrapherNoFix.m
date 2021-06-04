
Al=logspace(-4,1,100);
for j=1:100
    
I=eye(10);
U_Mat=zeros(10);
for i=1:10
[~,U_Mat(:,i),z]=StabCalc(I(:,i),0,CMatrix,SMatrix,Al(j));
%norm(U_Mat(:,i))/norm(I(:,i))
end
[u,s,v]=svd(U_Mat);
s=diag(s);
STAB(j)=max(s);
end
semilogx(Al,STAB)