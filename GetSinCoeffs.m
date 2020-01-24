
S0=FullSinCont*SinCoeffs(:,1);
S0=S0(8:10:end);
S1=FullSinCont*SinCoeffs(:,2);
S1=S1(8:10:end);
S2=FullSinCont*SinCoeffs(:,3);
S2=S2(8:10:end);
S3=FullSinCont*SinCoeffs(:,4);
S3=S3(8:10:end);
S4=FullSinCont*SinCoeffs(:,5);
S4=S4(8:10:end);
S5=FullSinCont*SinCoeffs(:,6);
S5=S5(8:10:end);
S6=FullSinCont*SinCoeffs(:,7);
S6=S6(8:10:end);
S7=FullSinCont*SinCoeffs(:,8);
S7=S7(8:10:end);
S8=FullSinCont*SinCoeffs(:,9);
S8=S8(8:10:end);
S9=FullSinCont*SinCoeffs(:,10);
S9=S9(8:10:end);

SMatrix=[S0,S1,S2,S3,S4,S5,S6,S7,S8,S9];

xx=linspace(-1,1,10);
[Q,R]=qr(fliplr(vander(xx)));
for i = 1:10
    Et(i)=max(abs(SMatrix(43:52,i)-Q(:,i)));
    if Et(i) > 1e-1
        Et(i)=max(abs(SMatrix(43:52,i)+Q(:,i)));
    end
end
SMatrix=[SMatrix(43:52,:);SMatrix(53:end,:);SMatrix(1:42,:)];

