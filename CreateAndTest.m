load AllMats.mat
n=length(CosCoeffs(1,:));
m=length(FullCosCont(1:9:end,1));

for i =1:n
    f=FullCosCont*CosCoeffs(:,i);
    f=f(1:9:end);
    ft=fft(f);
    ftm=[ft(1:m/2);zeros(9*m,1);ft(m/2+1:end)];
    ift=real(ifft(10*ftm));
    er(i)=norm(ift(1:10:end)-f);
end
