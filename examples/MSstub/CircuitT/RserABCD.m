function [T] = RserABCD(R)

% function [T] = RserABCD(R)
% Series resistor ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of R.
%
% Inputs:
% R - Resistance Ohm (can, in general, be function of frequency)

Nf = length(R);
T(1,1,:) = ones(Nf,1);
T(1,2,:) = R;
T(2,1,:) = zeros(Nf,1);
T(2,2,:) = T(1,1,:);

