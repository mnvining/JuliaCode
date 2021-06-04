clear all 
close all
GetCosPts
GetSinCoeffs
I=eye(10);
U_Mat=zeros(10);
for i=1:10
[~,U_Mat(:,i),z]=StabCalc(I(:,i),0,CMatrix,SMatrix,.001);
%norm(U_Mat(:,i))/norm(I(:,i))
end
[u,s,v]=svd(U_Mat);
s=diag(s);
bad=find(s>1);
n=length(bad);
VBad=v(:,bad);
%%
for i=1:10
    %[~,NewMat(:,i),z,R]=StabCalc(transpose(v(i,:)),1,CMatrix,SMatrix,.001);
    [~,U_Mat(:,i),z,R]=StabCalc(I(:,i),.14,CMatrix,SMatrix,.001);
    norm(U_Mat(:,i))/norm(I(:,i))
end

[u2,s2,v2]=svd(U_Mat);
[s,diag(s2)]
