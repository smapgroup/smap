function [Z0] = stripLine(eps_r,b,w,t)

% function [Z] = stripLine(eps_r,b,w,t)
% Calculates the characteristic impedance of symmetric stripline using the
% formulas in Balanis EM
% Much expansion possible!
% Inputs:
% eps_r - relative permittivity
% b - height of guide in m
% w - width of trace in m
% t - thickness of trace in m (can be zero)

load constants

tb = t./b;
tb1 = (1 - tb);
if t == 0
    Cf = eps0.*eps_r.*0.4413;
else
    Cf = eps0*eps_r./pi.*(2./tb1.*log(1 + 1./tb1) - (1./tb1 - 1).*log(1./tb1.^2 - 1));
end
Z0 = 1./sqrt(eps_r).*30.*pi./(w./(b*tb1) + Cf./(eps0*eps_r));