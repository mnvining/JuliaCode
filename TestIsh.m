
%N=round(N);
f=@(x)sqrt(1+x.^2);
thing=0:.1:10;
for i=1:length(thing)
    N=round(10*2^(thing(i)));
    x=[-1:2/(N-1):1]';
    x2=[-1:2/((N-1)*10):1]';
    y=f(x);
    plot(x,y)
    hold on
    ff=fcg(x,y,9,CMatrix,SMatrix);
    ft=fft(ff);
    Q=round((N+25)/2);
    pft=[ft(1:Q);zeros(9*(N+25),1);ft(Q+1:end)];
    ift=real(ifft(10*pft));
    nf=ift(1:10*(N-1)+1);
    er(i)=norm(nf-f(x2));
end

