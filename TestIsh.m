close all
%N=round(N);
f=@(x)exp(sin(2.7*pi*x)+cos(pi*x));
thing=0:.5:7;
NPlot=round(10.*2.^thing);
NHolder=[];
er=zeros(length(thing),1);
for i=1:length(thing)
    N=round(10*2^(thing(i)));
    NHolder=[NHolder,N];
    x=[-1:2/(N-1):1]';
    x2=[-1:2/((N-1)*10):1]';
    y=f(x);
    ff=fcg(y,5,CMatrix,SMatrix);
    xplot=[x;[1+2/(N-1):2/(N-1):1+50/(N-1)]'];
    ft=fft(ff);
    Q=round((N+25)/2);
    pft=[ft(1:Q);zeros(9*(N+25),1);ft(Q+1:end)];
    ift=real(ifft(10*pft));
    nf=ift(1:10*(N-1)+1);
    er(i)=max(abs((nf-f(x2))));
    
    if mod(i,5)==0
        figure
        plot(x,y,'b','Linewidth',3)
        hold on
        plot(xplot,ff,'k')
        title(['FC-Gram Continuation for f(x)=exp(sin(2.7pix)+cos(pix)), N=',num2str(N)])
    end
    
        
end
figure
loglog(NPlot,er)

