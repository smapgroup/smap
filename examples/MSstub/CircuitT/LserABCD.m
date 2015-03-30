function [T] = LserABCD(L,f)

% function [T] = LserABCD(L,f)
% Series inductance ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of f.
%
% Inputs:
% L - Inductance in H (can, in general, be function of frequency)
% f - frequency in Hz

Nf = length(f);
wL = 2.*pi.*f.*L;
T(1,1,:) = ones(Nf,1);
T(1,2,:) = 1i.*wL;
T(2,1,:) = zeros(Nf,1);
T(2,2,:) = T(1,1,:);

