function [Z] = ABCD2Z(ABCD)

% function [Z] = ABCD2S(ABCD)
% Converts the ABCD matrix to the Z matrix

A = ABCD(1,1,:);
B = ABCD(1,2,:);
C = ABCD(2,1,:);
D = ABCD(2,2,:);
Z(1,1,:) = A./C;
Z(1,2,:) = (A.*D - B.*C)./C;
Z(2,1,:) = 1./C;
Z(2,2,:) = D./C;