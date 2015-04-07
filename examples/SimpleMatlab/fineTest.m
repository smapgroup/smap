function Rf = fineTest(xi,f)
% f = linspace(0,pi,101);
x1 = 1.0*xi(1) + 0.0; 
x2 = 0.4*xi(2) - 0.0;
Rf = 1.1*sin(x1*f' + x2);
Rf = reshape(Rf,length(f),1);