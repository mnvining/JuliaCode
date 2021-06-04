function [u_sol,stab,acc]=SolCalc(fx,fs,N,c,CMatrix,SMatrix,Al)
x=[-1:2/(N-1):1]';
y=fx(x);
x2=[-1:2/((N-1)*10):1]';

f=fcgwc(y,9,c,CMatrix,SMatrix);
ft=fft(f);
Q=round((N+25)/2);
pft=[ft(1:Q);zeros(9*(N+25),1);ft(Q+1:end)];
Nx=length(pft);
k=2*pi*1i/((2*(N+25))/(N-1))*[[0:floor(Nx/2)]';[-floor(Nx/2)+1:-1]']; % 2*(N+25)/(N-1) is length of new interval with continuation
z=1-Al*k.^2;

AA=pft./z;
R=real(ifft(10*AA));


% orthogonalize solution to h1,h2 to check stability. if stability is
% achieved, all set.
h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));

[Q,~]=qr([h1,h2],0);
T=Q\R(1:10:10*N);
U=R(1:10:10*N)-Q*T;
stab=norm(U)/norm(y);

%create homogeneous on fine grid
h1f=exp((-x2-1)/sqrt(Al));
h2f=exp(-(1-x2)/sqrt(Al));

u_acc=fs(x2,Al);
A=[h1f(1) h2f(1); h1f(end) h2f(end)];
b=[R(1);R(10*(N-1)+1)];
c=A\b;
U_bc=R(1:10*(N-1)+1)-c(1)*h1f-c(2)*h2f;
%stab=norm(U)/norm(y);
acc=norm(U_bc-u_acc);
u_sol=U_bc(1:10:end);
