load('AllMats.mat')
C0=FullCosCont*CosCoeffs(:,1);
C0=C0(8:10:end);
C1=FullCosCont*CosCoeffs(:,2);
C1=C1(8:10:end);
C2=FullCosCont*CosCoeffs(:,3);
C2=C2(8:10:end);
C3=FullCosCont*CosCoeffs(:,4);
C3=C3(8:10:end);
C4=FullCosCont*CosCoeffs(:,5);
C4=C4(8:10:end);
C5=FullCosCont*CosCoeffs(:,6);
C5=C5(8:10:end);
C6=FullCosCont*CosCoeffs(:,7);
C6=C6(8:10:end);
C7=FullCosCont*CosCoeffs(:,8);
C7=C7(8:10:end);
C8=FullCosCont*CosCoeffs(:,9);
C8=C8(8:10:end);
C9=FullCosCont*CosCoeffs(:,10);
C9=C9(8:10:end);

CMatrix=[C0,C1,C2,C3,C4,C5,C6,C7,C8,C9];

xx=linspace(-1,1,10);
[Q,R]=qr(fliplr(vander(xx)));
for i = 1:10
    Et(i)=max(abs(CMatrix(43:52,i)-Q(:,i)));
    if Et(i) > 1e-1
        Et(i)=max(abs(CMatrix(43:52,i)+Q(:,i)));
    end
end
CMatrix=[CMatrix(43:52,:);CMatrix(53:end,:);CMatrix(1:42,:)];


