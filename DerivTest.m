N=10;
sc=10;
Al=1;
x=[-1:2/(N-1):1]';
y=x.^2-1;
f=fcg(x,y,9,CMatrix,SMatrix);
ft=fft(f);

pft=fft(f);
Nx=length(pft);
Q=round(Nx/2);

P=length(f);
k=2*pi/(70/9)*[[1:round(Nx/2)]';[-round(Nx/2)+1:-1]'];

z=1+Al*k.^2;

AA=pft./z;
R=real(ifft(AA));
plot(real(ifft(ft)))
hold on
plot(R)

stab=norm(R(1:10))/norm(y)
1-stab
h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));
A=[h1,h2];
b=R(1:10)-(y-2*Al);
c=A\b;
U=R(1:10)+c(1)*h1+c(2)*h2;
acc=abs(U-(y-2*Al))
hold off



