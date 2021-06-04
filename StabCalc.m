function [u_0c,u_qc,z,R]=StabCalc(fx,c,CMatrix,SMatrix,Al)
% [u_c]=StabCalc(fx,c,CMatrix,SMatrix,Al)
%               fx is vector.
N=length(fx);

f=fcgwc(fx,9,c,CMatrix,SMatrix);
ft=fft(f);
Q=round((N+25)/2);
pft=[ft(1:Q);zeros(9*(N+25),1);ft(Q+1:end)];
Nx=length(pft);
k=2*pi*1i/((2*(N+25))/(N-1))*[[0:floor(Nx/2)]';[-floor(Nx/2)+1:-1]']; % 2*(N+25)/(N-1) is length of new interval with continuation
z=1-Al*k.^2;

AA=pft./z;

R=real(ifft(10*AA));

x=[-1:2/(N-1):1]';
% orthogonalize solution to h1,h2 to check stability. if stability is
% achieved, all set.
h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));

[Q,~]=qr([h1,h2],0);
T=Q'*R(1:10:N*10);
U=R(1:10:N*10)-Q*T;
u_qc=U;
stab=norm(U)/norm(fx);


A=[h1(1) h2(1); h1(end) h2(end)];
b=[R(1);R(9*N+1)];
co=A\b;
u_0c=R(1:10:N*10)-co(1)*h1-co(2)*h2;
%stab=norm(U)/norm(y);
