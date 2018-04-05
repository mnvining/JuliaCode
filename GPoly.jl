function GPolyTest(Alpha,N,f,cl)
# function [n1,n2]=GPolyTest(Alpha,N,f,OPT)
# Alpha - alpha coeff of d2/dx2
# N - number of pts, 2^power please
# f - function handle for function to use
# cl - continuation length desired
# Out: n1 stability n2 accuracy

M=8
j=collect(2*[0:N-1]-N)
x=j./N
h=x[2]-x[1]

x2=collect(-1:h:1+cl-h)

u=f(x)
u=u'
fd=cl+cl/(N-1)

# Continuation Matrix:
Cont_Matrix=createFWDcont(N,vcat(collect(1:M+1),collect(2*N-M+1:2*N)),cl)
BackDouble=createBackDoubleDomain(N, vcat(collect(1:M+1),collect(2*N-M+1:2*N)),cl)

# Derivative Operator: (I-Alpha*d^2)
D=createId2Derivative(N,M,Alpha,cl)
D_Dag=inv(D)

# Compute solution G(conv)F to get Green's function solution
GG=GSol(Alpha,f,x);
[h1,h2,dH1dx,dH2dx]=CreateHSolOnly(Alpha,x);

[dGdx]=secder(GG,x);

Z=[dH1dx*(dH1dx') dH2dx*(dH1dx'); dH1dx*(dH2dx') dH2dx*(dH2dx')];

S=Z\[-dGdx*(dH1dx'); -dGdx*(dH2dx')];


Bf=GG+S(1)*h1+S(2)*h2;

figure(88)
plot(x,Bf,x,u,x,GG)
legend('Bf','f(x)','GG')


n1=norm(Bf);
n2=norm(u);

% Set-up matrix 
M=[Cont_Matrix;Cont_Matrix*D_Dag];

switch OPT
    case 0
        b=[u;GG'];
    case 1
        b=[u;Bf'];
end

% CIn is the c'*f that we need
CIn=M\b;



%% Norms and more Plotting
GreensCont=real(Cont_Matrix*CIn);
q=norm(GreensCont-u);
% q is the norm of the continuation vs 

sol=real(Cont_Matrix*D_Dag*CIn);
figure(44)
plot(x,sol)
Acc=norm(sol-(Bf'));
Stab=norm(sol)/norm(u);


solext=real(BackDouble*D_Dag*CIn);

figure(4)
plot(x2,solext,x2,real(BackDouble*CIn))
legend('deriv','orig')



