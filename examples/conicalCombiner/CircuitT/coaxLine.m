function [Z0] = coaxLine(r_in,r_out,eps_r)

% function [Z] = stripLine(r_in,r_out,eps_r)
% Calculates the characteristic impedance of coaxial line
% Much expansion possible!
% Inputs:
% r_in - innner conductor radius in m
% r_out - outer conductor radius in m
% eps_r - relative permittivity (optional)

load constants

Z0 = 60*log(r_out./r_in);