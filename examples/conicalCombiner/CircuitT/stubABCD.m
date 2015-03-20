function [T] = stubABCD(Z0,gamma,L,sp,os)

% function [T] = stubABCD(Z0,gamma,L,sp,os)
% Stub TXline ABCD matrix of size [2,2,Nf], where Nf is the number of 
% frequencies as reflected in the length of gamma.
%
% Inputs:
% Z0 - characteristic impedance in Ohm (can, in general, be function of frequency)
% gamma - complex propagation constant (is a real number is given, it will
%         be assumed that the line is lossless)
% L - Length in [m]
% sp - 's' for series and 'p' for shunt 
% os - 'o' for open and 's' for short

if isreal(gamma), gamma = 1i.*gamma; end % Check for "beta" propagation...

Nf = length(gamma);
gL = gamma.*L;

if isequal(os,'s')
    Zin = Z0.*tanh(gL);
elseif isequal(os,'o')
    Zin = Z0.*coth(gL);
end
if isequal(sp,'p')  % Shunt case
    T(1,2,:) = zeros(Nf,1);
    T(2,1,:) = 1./Zin;
elseif isequal(sp,'s')
    T(1,2,:) = Zin;
    T(2,1,:) = zeros(Nf,1);
end
% All cases
T(1,1,:) = ones(Nf,1);
T(2,2,:) = T(1,1,:);