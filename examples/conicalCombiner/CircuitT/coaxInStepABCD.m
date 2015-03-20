function [T] = coaxInStepABCD(di1,di2,do,f)

% function [T] = coaxInStepABCD(di1,di2,do,f)
% Coaxial line inner conductor step ABCD matrix of size [2,2,Nf], where Nf is the number of
% frequencies as reflected in the length of f.
%
% Inputs:
% di1 - Smaller inner diameter in m
% di2 - Larger inner diameter in m
% do - Outer diameter in m 
% f - frequency in Hz

load constants

if di1 > di2
    temp = di1;
    di1 = di2;
    di2 = temp;
end
alpha = (do - di2)/(do - di1);
tau = do/di1;

if alpha == 1
    C = 0;
else
    Cj = (eps0/(100*pi)*((alpha^2+1)/alpha*log((1+alpha)/(1-alpha)) - 2*log(4*alpha/(1-alpha^2))) + 1.11e-15*(1-alpha)*(tau-1));
    C = Cj*2*pi*(do/2)*100;
end
Nf = length(f);
wC = 2.*pi.*f.*C;
T(1,1,:) = ones(Nf,1);
T(1,2,:) = zeros(Nf,1);
T(2,1,:) = 1i.*wC;
T(2,2,:) = T(1,1,:);

