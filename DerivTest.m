close all
N=10;
sc=10;
Al=1;
x=[-1:2/(N-1):1]';
y=x.^2;
f=fcg(y,9,CMatrix,SMatrix);
ft=fft(f);
yhat=fft(y);

pft=fft(f);
Nx=length(pft);
Q=round(Nx/2);

P=length(f);
k=2*pi*1i/(70/9)*[[0:floor(Nx/2)]';0;[-floor(Nx/2)+1:-1]'];
td=k.^2.*ft;
figure
i2=real(ifft(td));
plot(i2)

 
z=1-Al*k.^2;

AA=pft./z;
R=real(ifft(AA));

stab=norm(R(1:10))/norm(y)
norm(y+2*Al)/norm(y)
1-stab
h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));
A=[h1(1) h2(1); h1(end) h2(end)]
b=[R(1);R(10)]
c=A\b;

U=R(1:10)-c(1)*h1-c(2)*h2;
plot(U)
norm(U)/norm(y)


