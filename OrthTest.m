c=0;
Al=.0001;
N=100;
x=[-1:2/(N-1):1]';
%fx=@(x)cos(100/7*x)
%fx=@(x)sin(2*x);
fx=@(x)exp(x);
%fx=@(x)x.^2;
%fx=@(x)1./(1+25*x.^2);
%y=x.^2;
y=exp(x);
%y=sin(2*x);
%y=cos(100/7*x);
% runge function <-- good test
% e^(cosx) or e^10x<- BE CAREFUL W ALPHA
%y=1./(1+25*x.^2);


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
plot(R)


a=Al;
x_c=-1:2/(N-1):1;
h1=exp((-x-1)/sqrt(Al));
h2=exp(-(1-x)/sqrt(Al));
%c1=-(1+2*Al)*exp(2/sqrt(Al))/(exp(2/sqrt(Al))+1); % note c1=c2;
%u_acc=x.^2+2*Al+c1*h1+c1*h2;
% for runge: u_acc= -(1i*exp(-x/sqrt(a)-1i/(5*sqrt(a))).*(expint(-(1 + 1i/5)/sqrt(a))*exp((2*x)/sqrt(a)+(2 + (2*1i)/5)/sqrt(a))-expint(-(1 - 1i/5)/sqrt(a))*exp((2*x)/sqrt(a)+2/sqrt(a))+expint((1 - 1i/5)/sqrt(a))*exp((2*x)/sqrt(a) + (2*1i)/(5*sqrt(a)))-expint((1 + 1i/5)/sqrt(a))*exp((2*x)/sqrt(a))-exp((2*x)/sqrt(a) + (2*1i)/(5*sqrt(a)))*expint((-5*x - 1i)/(5*sqrt(a)))- exp((2*x)/sqrt(a) + (2 + (2*1i)/5)/sqrt(a))*expint((-5*x - 1i)/(5*sqrt(a))) + exp((2*x)/sqrt(a))*expint((1i - 5*x)/(5*sqrt(a))) + exp((2*x)/sqrt(a) + 2/sqrt(a))*expint((1i - 5*x)/(5*sqrt(a))) - exp((2*1i)/(5*sqrt(a)))*expint((5*x - 1i)/(5*sqrt(a))) - exp((2 + (2*1i)/5)/sqrt(a))*expint((5*x - 1i)/(5*sqrt(a))) + exp(2/sqrt(a)).* expint((5*x + 1i)/(5*sqrt(a))) + expint((5*x + 1i)/(5*sqrt(a))) + exp((2 + (2*1i)/5)/sqrt(a))*expint(-(1 + 1i/5)/sqrt(a)) - exp(2/sqrt(a)).* expint(-(1 - 1i/5)/sqrt(a)) + exp((2*1i)/(5*sqrt(a))).* expint((1 - 1i/5)/sqrt(a)) - expint((1 + 1i/5)/sqrt(a))))/(20*sqrt(a).*(exp(2/sqrt(a)) + 1));
%u_acc= -(1i*exp(-x/sqrt(a)-1i/(5*sqrt(a))).*(expint(-(1 + 1i/5)/sqrt(a)).*exp((2*x)/sqrt(a)+(2 + (2*1i)/5)/sqrt(a))-expint(-(1 - 1i/5)/sqrt(a)).*exp((2*x)/sqrt(a)+2/sqrt(a))+expint((1 - 1i/5)/sqrt(a)).*exp((2*x)/sqrt(a) + (2*1i)/(5*sqrt(a)))-expint((1 + 1i/5)/sqrt(a)).*exp((2*x)/sqrt(a))-exp((2*x)/sqrt(a) + (2*1i)/(5*sqrt(a))).*expint((-5*x - 1i)/(5*sqrt(a)))- exp((2*x)/sqrt(a) + (2 + (2*1i)/5)/sqrt(a)).*expint((-5*x - 1i)/(5*sqrt(a))) + exp((2*x)/sqrt(a)).*expint((1i - 5*x)/(5*sqrt(a))) + exp((2*x)/sqrt(a) + 2/sqrt(a)).*expint((1i - 5*x)/(5*sqrt(a))) - exp((2*1i)/(5*sqrt(a))).*expint((5*x - 1i)/(5*sqrt(a))) - exp((2 + (2*1i)/5)/sqrt(a)).*expint((5*x - 1i)/(5*sqrt(a))) + exp(2/sqrt(a)).* expint((5*x + 1i)/(5*sqrt(a))) + expint((5*x + 1i)/(5*sqrt(a))) + exp((2 + (2*1i)/5)/sqrt(a)).*expint(-(1 + 1i/5)/sqrt(a)) - exp(2/sqrt(a)).* expint(-(1 - 1i/5)/sqrt(a)) + exp((2*1i)/(5*sqrt(a))).* expint((1 - 1i/5)/sqrt(a)) - expint((1 + 1i/5)/sqrt(a))))./(20*sqrt(a).*(exp(2/sqrt(a)) + 1));
%u_acc=-exp(x2/sqrt(Al))*exp(-2/(5*sqrt(-Al)))*exp(-(3*sqrt(-Al))/(5*Al))*expint(1, -1/sqrt(Al) + 1/(5*sqrt(-Al)))/(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) - exp(x2/sqrt(Al))*exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al))^2*exp(-(3*sqrt(-Al))/(5*Al))*expint(1, 1/sqrt(Al) + 1/(5*sqrt(-Al)))/(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(x2/sqrt(Al))*exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al))^2*exp(-sqrt(-Al)/(5*Al))*expint(1, 1/sqrt(Al) - 1/(5*sqrt(-Al)))/(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(x2/sqrt(Al))*exp(-2/(5*sqrt(-Al)))*expint(1, -1/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(-sqrt(-Al)/(5*Al))/(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) - exp(-2/(5*sqrt(-Al)))*exp(-(3*sqrt(-Al))/(5*Al))*expint(1, -1/sqrt(Al) + 1/(5*sqrt(-Al)))/(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al))-exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al))^2*exp(-(3*sqrt(-Al))/(5*Al))*expint(1, 1/sqrt(Al) + 1/(5*sqrt(-Al)))/(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al))^2*exp(-sqrt(-Al)/(5*Al))*expint(1, 1/sqrt(Al) - 1/(5*sqrt(-Al)))/(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al)))*expint(1, -1/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(-sqrt(-Al)/(5*Al))/(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al)))*expint(1, -x2/sqrt(Al) + 1/(5*sqrt(-Al)))*exp(-(3*sqrt(-Al))/(5*Al))/(20*sqrt(-Al)*exp(x2/sqrt(Al))) + exp(-2/(5*sqrt(-Al)))*expint(1, x2/sqrt(Al) + 1/(5*sqrt(-Al)))*exp(-(3*sqrt(-Al))/(5*Al))*exp(x2/sqrt(Al))/(20*sqrt(-Al)) - exp(-2/(5*sqrt(-Al)))*expint(1, -x2/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(-sqrt(-Al)/(5*Al))/(20*sqrt(-Al)*exp(x2/sqrt(Al)))-exp(-2/(5*sqrt(-Al)))*expint(1, x2/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(x2/sqrt(Al))*exp(-sqrt(-Al)/(5*Al))/(20*sqrt(-Al));

% prior to homogeneous work, orthogonalize h1,h2 1st to each other [QR]
% then project into soln to orthogonalize soln to u_h.
[Q,~]=qr([h1,h2],0);
%Q1=Q(:,1); % q1 is orthogonalized h1/h2 basis1
%Q2=Q(:,2); % q2 is orthogonalized h1/h2 basis2

T=Q\R(1:10:10*N);

U=R(1:10:10*N)-Q*T;

stab=norm(U)/norm(y)




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





h1f=exp((-x2-1)/sqrt(Al));
h2f=exp(-(1-x2)/sqrt(Al));
%c1=-(1+2*Al)*exp(2/sqrt(Al))/(exp(2/sqrt(Al))+1); % note c1=c2;
%u_acc=fx(x2)+2*Al+c1*h1f+c1*h2f;
% c1=(exp(1-2/sqrt(Al))-exp(-1))/((Al-1)(exp(-4/sqrt(Al))-1)
%c1 = (exp((sqrt(Al) - 2)/sqrt(Al)) - exp(-1))/((-1 + Al)*(exp(-4/sqrt(Al)) - 1));
% c2=(exp(-1-2/sqrt(Al))+exp(1))/((Al-1)(exp(-4/sqrt(Al))-1)
%c2=(exp(-(2 + sqrt(Al))/sqrt(Al)) - exp(1))/((-1 + Al)*(exp(-4/sqrt(Al)) - 1));
%u_acc=1/(1-Al)*fx(x2)+c1*h1f+c2*h2f;
u_acc=(exp(-(x2 + sqrt(Al) - 1)./sqrt(Al)) - exp((-x2 - 1 + sqrt(Al))./sqrt(Al)) - exp((x2 - sqrt(Al) - 1)./sqrt(Al)) + exp((x2 + 1 + sqrt(Al))/sqrt(Al)) - exp((x2*sqrt(Al) + 2)./sqrt(Al)) + exp((x2*sqrt(Al) - 2)./sqrt(Al))).*exp(2/sqrt(Al))./((Al - 1)*(exp(4/sqrt(Al)) - 1));
%u_acc=(exp((2 + x2)./sqrt(Al)).*sin(2*x2) - exp((2*x2 + 1)./sqrt(Al)).*sin(2) + exp(1/sqrt(Al))*sin(2) - exp(x2./sqrt(Al)).*sin(2*x2)).*exp(-x2/sqrt(Al))./((exp(2/sqrt(Al)) - 1)*(1 + 4*Al));
%u_acc=(exp((5*sqrt(Al) + 1 + x2)./sqrt(Al)) - exp(-(-5*sqrt(Al) + 1 + x2)./sqrt(Al)) + exp(-(5*sqrt(Al) - 1 + x2)./sqrt(Al)) - exp((-5*sqrt(Al) - 1 + x2)./sqrt(Al)) - exp((5*x2*sqrt(Al) + 2)/sqrt(Al)) + exp((5*x2*sqrt(Al) - 2)/sqrt(Al))).*exp(2/sqrt(Al))/((25*Al - 1).*(exp(4/sqrt(Al)) - 1));
%u_acc=-49*(exp((2*x2 + 1)/sqrt(Al)).*cos(100/7) - exp((2 + x2)/sqrt(Al)).*cos((100*x2)/7) + exp(1/sqrt(Al))*cos(100/7) - exp(x2/sqrt(Al)).*cos((100*x2)/7)).*exp(-x2/sqrt(Al))./((exp(2/sqrt(Al)) + 1)*(49 + 10000*Al));
%u_acc=-exp(x2/sqrt(Al)).*exp(-2/(5*sqrt(-Al))).*exp(-(3*sqrt(-Al))./(5*Al)).*expint(-1/sqrt(Al) + 1/(5*sqrt(-Al)))./(20*(exp(1/sqrt(Al)).^2 + 1)*sqrt(-Al)) - exp(x2/sqrt(Al)).*exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al)).^2.*exp(-(3*sqrt(-Al))/(5*Al)).*expint(1/sqrt(Al) + 1/(5*sqrt(-Al)))./(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(x2/sqrt(Al)).*exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al)).^2.*exp(-sqrt(-Al)/(5*Al)).*expint(1/sqrt(Al) - 1/(5*sqrt(-Al)))./(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(x2/sqrt(Al)).*exp(-2/(5*sqrt(-Al))).*expint(-1/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(-sqrt(-Al)/(5*Al))/(20*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) - exp(-2/(5*sqrt(-Al)))*exp(-(3*sqrt(-Al))/(5*Al)).*expint(-1/sqrt(Al) + 1/(5*sqrt(-Al)))./(20.*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al)).^2 + 1)*sqrt(-Al))-exp(-2/(5*sqrt(-Al))).*exp(1/sqrt(Al)).^2*exp(-(3*sqrt(-Al))/(5*Al)).*expint(1/sqrt(Al) + 1/(5*sqrt(-Al)))./(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al)))*exp(1/sqrt(Al)).^2.*exp(-sqrt(-Al)/(5*Al)).*expint(1/sqrt(Al) - 1/(5*sqrt(-Al)))./(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al))).*expint(-1/sqrt(Al) - 1/(5*sqrt(-Al))).*exp(-sqrt(-Al)/(5*Al))/(20*exp(x2/sqrt(Al)).*(exp(1/sqrt(Al))^2 + 1)*sqrt(-Al)) + exp(-2/(5*sqrt(-Al))).*expint(-x2/sqrt(Al) + 1/(5*sqrt(-Al))).*exp(-(3*sqrt(-Al))/(5*Al))./(20*sqrt(-Al).*exp(x2/sqrt(Al))) + exp(-2/(5*sqrt(-Al))).*expint(x2/sqrt(Al) + 1/(5*sqrt(-Al))).*exp(-(3*sqrt(-Al))/(5*Al)).*exp(x2/sqrt(Al))./(20*sqrt(-Al)) - exp(-2/(5*sqrt(-Al))).*expint(-x2/sqrt(Al) - 1/(5*sqrt(-Al)))*exp(-sqrt(-Al)./(5*Al))./(20*sqrt(-Al).*exp(x2/sqrt(Al)))-exp(-2/(5*sqrt(-Al))).*expint(x2/sqrt(Al) - 1/(5*sqrt(-Al))).*exp(x2/sqrt(Al)).*exp(-sqrt(-Al)/(5*Al))./(20*sqrt(-Al));


A=[h1f(1) h2f(1); h1f(end) h2f(end)];
b=[R(1);R(10*(N-1)+1)];
c=A\b;
U_bc=R(1:10*(N-1)+1)-c(1)*h1f-c(2)*h2f;
%stab=norm(U)/norm(y);
norm(U_bc-u_acc)

% temp assign uc, uf
u_c=U;
u_f=0;


% function inputs: 
% stab parameter
% return the coarse solution
% u_c, u_f

