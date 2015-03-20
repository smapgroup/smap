function x = decReal(c,Mvect,xmin,xmax)
% function x = decReal(c,Mvect,xmin,xmax)
% decodes the chromosome c to the real vector x.
% numel(Mvect) = numel(xmin) = numel(xmax) = numel(x)

Nx = numel(Mvect);   % Number of variables
if numel(xmin) ~= numel(xmax)
    error('xmax should be the same length as xmin');
elseif Nx ~= numel(xmin)
    error('Mvect should be the same length as xmin and xmax');
end

for ii = 1:Nx
    % Get the current variable chromosome
    if ii == 1
        ci = c(1:Mvect(1));
    else
        ci = c(sum(Mvect(1:ii-1))+1:sum(Mvect(1:ii)));
    end
    normBin = 2.^Mvect(ii) - 1;
    ciStr = num2str(ci);
    ciS = strrep(ciStr,' ',''); % Remove white spaces
    xpu = bin2dec(ciS)./normBin;
    x(ii) = (xmax(ii) - xmin(ii))*xpu + xmin(ii);
end
