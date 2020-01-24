function C=ptc(y,d,CMatrix)
% function C=ptstocoeffs(x,y)
% calculates the coefficients of the gram polynomial expansion for the
% points (x,y).
% inputs: x - x-coords
%           y - corresponding y-coords [must be column vector]
% outputs: c - coefficients c_0,c_1... corresponding to the Gram Polys.

% Uses 0 - d.


% loads the even coefficients for degrees 0:8
p=CMatrix(1:10,1:d+1);
size(p)
size(y)
C=p\y;