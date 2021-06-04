function [u_c,u_f,stab]=TBD(y,Al,bl,br,xl,xr,c,CMatrix,SMatrix)
% function function [u_c,u_f]=TBD(y,Al,bl,br,xl,xr,c,CMatrix,SMatrix)
% Inputs:       y - function to be approximated at N points
%               Al - alpha, parameter in (u-alphau_xx)=f
%               bl,br - boundary values on the left, and right
%               xl,xr - whether or not to start at the boundary or on the
%               boundary
%               c - continuation scaling parameter
%               CMatrix/SMatrix - coeffs for even/odd continuations 
% Outputs:      u_c - solution u evaluated on coarse grid
%               u_f - solution u evaluated on the fine grid 

N=length(y);



f=fcg(y,9,CMatrix,SMatrix);
ft=fft(f);
Nx=length(ft);
k=2*pi*1i/(70/9)*[[0:floor(Nx/2)]';0;[-floor(Nx/2)+1:-1]'];
z=1-Al*k.^2;

AA=ft./z;
R=real(ifft(AA));



h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));

% prior to homogeneous work, orthogonalize h1,h2 1st to each other [QR]
% then project into soln to orthogonalize soln to u_h.
[Q,~]=qr([h1',h2'],0);
%Q(:,1); % q1 is orthogonalized h1/h2 basis1
%Q(:,2); % q2 is orthogonalized h1/h2 basis2

T=Q\R(1:10);

U=R(1:10)-Q*T;
plot(U)
% determine if the solution can be stabilized by orthogonalizing wrt homogeneous solutions

norm(U)/norm(y) 


% In thesis: prove that stab on 0 is trivial 

% Systematic study: check stability with c:
%                   what is parameter for 
%                                       line, quadratic, etc for 10 pts
%                   just discrete 10,20,100 pts
%                   take orth basis on ~25 pts (QR)
%                               Act on ID one at a time, generate operator
%                               acting on column to get the matrix
%                               representation (dep on alpha). 
%                                       check SVD to get norm
%                   check exp(x) on coarse grid                          




if xl==0
    Lind=1;
else
    Lind=2;
end
if xr==0
    Rind=N+1;
else
    Rind=N;
end




A=[h1(Lind) h2(Lind); h1(Rind) h2(Rind)];
b=[R(Lind)+bl;R(Rind)+br];
c=A\b;
U=R(1:10)-c(1)*h1-c(2)*h2;
stab=norm(U)/norm(y);

% temp assign uc, uf
u_c=U;
u_f=0;


% function inputs: 
% stab parameter
% return the coarse solution
% u_c, u_f

