function f=fcgram(x,y,d,CMatrix,SMatrix)
% function f=fcgram(x,y)
% input:    x, on [-1,1]
%           y, on [-1,1]
% First 8 y-values
y1=y(1:10);
% Last 8 y-values
y2=y(end-9:end);

% on the right side
C1=ptc(y2,d,CMatrix);
% on the left side
C2=ptc(y1,d,CMatrix);


% loads the even coefficients for degrees 0:9
% CMatrix
% loads the odd coefficients for degrees 0:9
% SMatrix


LSE=CMatrix(:,1:d+1)*C1;
LSO=SMatrix(:,1:d+1)*C1;
RSE=CMatrix(:,1:d+1)*C2;
RSO=SMatrix(:,1:d+1)*C2;

f=RSE/2-RSO/2+(LSE/2+LSO/2);


