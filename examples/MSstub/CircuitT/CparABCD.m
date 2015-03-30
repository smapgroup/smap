function [T] = CparABCD(C,f)

% function [T] = CparABCD(L,f)
% Shunt capacitance ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of f.
%
% Inputs:
% C - Capacitence in F (can, in general, be function of frequency)
% f - frequency in Hz

Nf = length(f);
wC = 2.*pi.*f.*C;
T(1,1,:) = ones(Nf,1);
T(1,2,:) = zeros(Nf,1);
T(2,1,:) = 1i.*wC;
T(2,2,:) = T(1,1,:);

