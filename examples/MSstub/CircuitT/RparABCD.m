function [T] = RparABCD(R)

% function [T] = RparABCD(R)
% Shunt resistor ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of R.
%
% Inputs:
% R - Resistance Ohm (can, in general, be function of frequency)

Nf = length(R);
T(1,1,:) = ones(Nf,1);
T(1,2,:) = zeros(Nf,1);
T(2,1,:) = 1./R;
T(2,2,:) = T(1,1,:);

