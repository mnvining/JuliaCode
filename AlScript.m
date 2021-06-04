
fx1=@(x)cos(100/7*x);
fx2=@(x)sin(2*x);
fx3=@(x)exp(x);
fx4=@(x)x.^2;

u_acc1=@(x2,Al)-49*(exp((2*x2 + 1)/sqrt(Al)).*cos(100/7) - exp((2 + x2)/sqrt(Al)).*cos((100*x2)/7) + exp(1/sqrt(Al))*cos(100/7) - exp(x2/sqrt(Al)).*cos((100*x2)/7)).*exp(-x2/sqrt(Al))./((exp(2/sqrt(Al)) + 1)*(49 + 10000*Al));
u_acc2=@(x2,Al)(exp((2 + x2)./sqrt(Al)).*sin(2*x2) - exp((2*x2 + 1)./sqrt(Al)).*sin(2) + exp(1/sqrt(Al))*sin(2) - exp(x2./sqrt(Al)).*sin(2*x2)).*exp(-x2/sqrt(Al))./((exp(2/sqrt(Al)) - 1)*(1 + 4*Al));
u_acc3=@(x2,Al)(exp(-(x2 + sqrt(Al) - 1)./sqrt(Al)) - exp((-x2 - 1 + sqrt(Al))./sqrt(Al)) - exp((x2 - sqrt(Al) - 1)./sqrt(Al)) + exp((x2 + 1 + sqrt(Al))/sqrt(Al)) - exp((x2*sqrt(Al) + 2)./sqrt(Al)) + exp((x2*sqrt(Al) - 2)./sqrt(Al))).*exp(2/sqrt(Al))./((Al - 1)*(exp(4/sqrt(Al)) - 1));
c1=@(Al)-(1+2*Al)*exp(2/sqrt(Al))/(exp(2/sqrt(Al))+1); % note c1=c2;
h1f=@(x2,Al)exp((-x2-1)/sqrt(Al));
h2f=@(x2,Al)exp(-(1-x2)/sqrt(Al));
u_acc4=@(x2,Al)fx4(x2)+2*Al+c1(Al)*h1f(x2,Al)+c1(Al)*h2f(x2,Al);


Al=.1;
N=2*round(10.^(linspace(1,4,100)));
N(1:5)
s=zeros(length(N),4);
a=s;
for i=1:length(N)
    [~,s(i,1),a(i,1)]=SolCalc(fx1,u_acc1,N(i),0,CMatrix,SMatrix,Al);
    [~,s(i,2),a(i,2)]=SolCalc(fx2,u_acc2,N(i),0,CMatrix,SMatrix,Al);
    [~,s(i,3),a(i,3)]=SolCalc(fx3,u_acc3,N(i),0,CMatrix,SMatrix,Al);
    [~,s(i,4),a(i,4)]=SolCalc(fx4,u_acc4,N(i),0,CMatrix,SMatrix,Al);
end
loglog(N,a)
axis tight
legend('y=cos(100/7x','y=sin(2x)','y=exp(x)','y=x^2','Location','Best') 
title('Accuracy based on number of points N, \alpha=0.1')
xlabel('N')
ylabel('norm(u_c - exact)')

% look @ convergence as a fn of N; 
% SVD of a matrix use V as the inputs; 
% use I
% show singular vectors can be stabilized



