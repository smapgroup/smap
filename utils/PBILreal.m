function [x,fval,exitFlag,output] = PBILreal(funfcn,xmin,xmax,M,options)

% function [x,fval,exitFlag,output] = PBIL(funfcn,xmin,xmax,M,options)
% Wrapper for PBIL which accepts an array of real input variables - this
% one is used most often of all the wrappers for PBIL...
%
% Returns the same as PBIL, except the cromosome c is replaced by the real
% values x.
% Output has an extra entry: 
% output.XFbest = XFbest;    % Best position and fuction value at each iteration [1,Nx+1]
%
% Input vectors xmin and xmax are row vectors containing the minimum and
% maximum values af the parameter ranges
% Input vector M is the number of bits to use for each variable.  If M is a
% scalar the same value will be used for all varables.  Default 8.
%
% See PBIL.m for more info...
% 
% Date created: 2014-11-13
% Dirk de Villiers
% Last Modified: 2014-11-13
% Updates:
% 2014-11-13: Build function

if nargin < 4
    M = 8;
    options = [];
elseif nargin < 5
    options = [];
end

if length(xmin) ~= length(xmax)
    error('xmax should be the same length as xmin');
else
    Nx = length(xmin);
end

if Nx > 1 && length(M) ~= Nx
    M = repmat(M(1),1,Nx);
end
xmin = reshape(xmin,1,Nx);
xmax = reshape(xmax,1,Nx);
% [c,fval,exitFlag,output] = PBIL(wrapper,sum(M),options);
[c,fval,exitFlag,output] = PBIL(@(c) wrapper(c,funfcn,M,xmin,xmax),sum(M),options);

output.XFbest = zeros(output.iterCount,Nx+1);
for ii = 1:output.iterCount
   output.XFbest(ii,1:Nx) = decReal(output.CFbest(ii,1:sum(M)),M,xmin,xmax);
   output.XFbest(ii,Nx+1) = output.CFbest(ii,sum(M)+1);
end
x = output.XFbest(end,1:Nx);
end


function fx = wrapper(c,realFunc,Mvect,xmin,xmax)
x = decReal(c,Mvect,xmin,xmax);
fx = realFunc(x);
end
