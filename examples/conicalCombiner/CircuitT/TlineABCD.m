function [T] = TlineABCD(Z0,gamma,L)

% function [T] = TlineABCD(Z0,gamma,L)
% Transmission line ABCD matrix of size [2,2,Nf], where Nf is the number of 
% frequencies as reflected in the length of gamma.
%
% Inputs:
% Z0 - characteristic impedance in Ohm (can, in general, be function of frequency)
% gamma - complex propagation constant (is a real number is given, it will
%         be assumed that the line is lossless)
% L - Length in [m]

if isreal(gamma), gamma = 1i.*gamma; end % Check for "beta" propagation...

gL = gamma.*L;
T(1,1,:) = cosh(gL);
T(1,2,:) = Z0.*sinh(gL);
T(2,1,:) = 1./Z0.*sinh(gL);
T(2,2,:) = T(1,1,:);

