using GenericSVD
i=5
N=2^(i+2)
K=2*N;
M=Int(N/4)
f = BigFloat(2*pi/K);
w = (0:f:2*pi-f/2)' * im;
x = 0:K-1;
qq=-w.*x;
e=exp(1)
A = 2^(i-4)/M*e.^((qq));
C=A[1:N,[1:M+1,end-M+1:end]];
Cinv=A[1:end,[1:M+1;end-M+1:end]];
(Uc,Sc,Vc)=svd(C);
Sc=convert(Array{AbstractFloat},Sc)
KK=-pi*im*[0:N-1;0;-(N-1):-1]';
n=vcat(KK[1:M+1],KK[(end-M+1):end]);
b=BigFloat(0.0001)
M=eye(length(n))-diagm(b*n.*n);
Sps=Sc;
Sps[Sps[:].>0]=1./Sps[Sps[:].>0]
Mp=convert(Array{AbstractFloat},M);
Mp[Mp.>0]=1./Mp[Mp.>0]
Op=C*Mp*Vc*diagm(Sps)'*Uc'
(Uo,So,Vo)=svd(Op)
for j=1:length(So)
    println(So[j])
end
numg1=sum(So.>1)





function 
