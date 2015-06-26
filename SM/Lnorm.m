function Ln = Lnorm(y,L)
% Calculates the norm of vector y.  L = 'L1' returns the L1 norm and L
% = 'L2' the L2 norm
if strcmp(L,'L1'), Lp = 1;
elseif strcmp(L,'L2'), Lp = 2; end
Ln = sum(abs(y).^Lp)./length(y);
end