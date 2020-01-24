N=logspace(1,4,50);

for i=1:50
    x=linspace(-1,1,N(i));
    y=exp(-x');
    ff=fcg(x,y,9,CMatrix,SMatrix);
    figure
    plot(ff)
    hold on
    plot(y)
end