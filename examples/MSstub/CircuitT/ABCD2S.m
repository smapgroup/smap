function [S] = ABCD2S(ABCD,Z01,Z02)

% function [S] = ABCD2S(ABCD,Z01,Z02)
% Converts the ABCD matrix to the S matrix in impedance Z01 (port 1) and Z02 (port 2) 
% Z1/2 can, in general, be a function of frequency with size Nf the same as
% the third dimension of ABCD [2,2,Nf]

if nargin < 3, Z02 = Z01; end

A = ABCD(1,1,:);
B = ABCD(1,2,:);
C = ABCD(2,1,:);
D = ABCD(2,2,:);

den = A.*Z02 + B + C.*Z01.*Z02 + D.*Z01;
S11 = (A.*Z02 + B - C.*conj(Z01).*Z02 - D.*conj(Z01))./den;
S12 = (2.*(A.*D - B.*C).*sqrt(real(Z01).*real(Z02)))./den;
S21 = (2.*sqrt(real(Z01).*real(Z02)))./den;
S22 = (-A.*conj(Z02) + B - C.*Z01.*conj(Z02) + D.*Z01)./den;
S(1,1,:) = S11;
S(1,2,:) = S12;
S(2,1,:) = S21;
S(2,2,:) = S22;