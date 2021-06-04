function f=fcgwc(y,d,c,CMatrix,SMatrix)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
% First 8 y-values
y1=y(10:-1:1);
% Last 8 y-values
y2=y(end-9:1:end);

% on the right side
C1=ptc(y2,d,CMatrix).*((1-c).^[0:9])';
% on the left side
C2=ptc(y1,d,CMatrix).*((1-c).^[0:9])';


% loads the even coefficients for degrees 0:9
% CMatrix
% loads the odd coefficients for degrees 0:9
% SMatrix

E1=(C1+C2)/2;
O1=(C1-C2)/2;


LSE=CMatrix(11:35,1:d+1)*E1;
LSO=SMatrix(11:35,1:d+1)*O1;
%RSE=CMatrix(11:35,1:d+1)*0;
%RSO=SMatrix(11:35,1:d+1)*0;

%f=-(RSE/2-RSO/2-(LSE/2+LSO/2));
% w=[1:25]';
% plot(w,LSE,w,LSO)
f=[y;LSE+LSO];

