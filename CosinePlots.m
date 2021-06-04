GetCosPts
GetSinCoeffs
h=2/9;
xPlot=-1:h:70*h-1-h;
for i=1:10
    figure(i)
    plot(xx,SMatrix(1:10,i),'r','linewidth',3)
    hold on
    plot(xPlot,SMatrix(:,i))
    plot(xPlot(36:45),SMatrix(36:45,i),'r','linewidth',3)
    xticks([-1:1:14])
    axis tight
    legend('Original Function','Continuation')
    title([num2str(i-1),'^{th} Degree Polynomial and Sine Continuation'])
end

    