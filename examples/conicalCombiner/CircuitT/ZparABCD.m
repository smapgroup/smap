function [T] = ZparABCD(Z)

% function [T] = ZparABCD(Z)
% Shunt impedance ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of Z.
%
% Inputs:
% Z - Impedance Ohm (can, in general, be function of frequency)

Nf = length(Z);
T(1,1,:) = ones(Nf,1);
T(1,2,:) = zeros(Nf,1);
T(2,1,:) = 1./Z;
T(2,2,:) = T(1,1,:);

