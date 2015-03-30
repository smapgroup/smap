function [Zf,beta] = ustripLine(eps_r,h,w,t,f)

% function [Z] = ustripLine(eps_r,h,w,t,f)
% Calculates the characteristic impedance and propagation constant
% of microstrip line using the formulas in Balanis EM
% Much expansion possible!
% Inputs:
% eps_r - relative permittivity
% h - height of line in m
% w - width of trace in m
% t - thickness of trace in m (can be zero)
% f - frequency in Hz (can be vector)

load constants

w_h = w./h;
t_h = t./h;
w_eff = w_h;
if w_h >= 1/(2*pi)
    if t_h > 0
        w_eff = w_eff + 1.25./pi.*t_h.*(1 + log(2./t_h));
    end
else
    if t_h > 0
        w_eff = w_eff + 1.25./pi.*t_h.*(1 + log(4*pi.*w./t));
    end
end
w_eff = w_eff.*h;

if w_eff/h <= 1
    eps_r_eff = (eps_r + 1)/2 + (eps_r - 1)./2.*((1 + 12.*h./w_eff).^-0.5 + 0.04.*(1 - w_eff./h).^2);
    Z0 = 60./sqrt(eps_r_eff).*log(8*h./w_eff + w_eff./(4.*h));
else
    eps_r_eff = (eps_r + 1)/2 + (eps_r - 1)./2.*(1 + 12.*h./w_eff).^-0.5;
    Z0 = 120.*pi./sqrt(eps_r_eff)./(w_eff./h + 1.393 + 0.667.*log(w_eff./h + 1.444));
end
    
lam0 = c0./f;
ft = Z0./(2.*mu0.*h);
eps_r_eff_f = eps_r - (eps_r - eps_r_eff)./(1 + eps_r_eff./eps_r.*(f./ft).^2);
lam_g = lam0./sqrt(1.*eps_r_eff_f);
beta = 2.*pi./lam_g;
Zf = Z0.*sqrt(eps_r_eff./eps_r_eff_f);
