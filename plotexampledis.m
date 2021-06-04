close all
f=@(x)1./(1+12*x.^2);
N=10;

x=[-1:2/(N-1):1]';
x2=[-1:2/(N-1):1+50/(N-1)]';
y=f(x);

[f,le,lo,re,ro]=fcgwplotter(y,9,CMatrix,SMatrix);
re=re(1:end-10);
le=le(11:end);
lo=lo(11:end);
ro=ro(1:end-10);
figure
plot([y;re])
hold on
plot([y;ro])
figure
plot([y;le])
hold on
plot([y;lo])
plot([y;(le+lo)])
figure

plot(x2,[y;(le+lo)/2-(re-ro)/2])


% figure
% plot(x,y,'b','Linewidth',4)
% hold on
% plot(x2,(le+re(end:-1:1))/2)
% legend('Original Function, last n points','Even continuation')
% 
% figure
% plot(x,y,'b','Linewidth',4)
% hold on
% plot(x2,(lo+ro(end:-1:1))/2)
% legend('Original Function, first n points (backwards)','Odd continuation')