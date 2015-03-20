function R = conicalCoarseModelDdV(x,xp,f)
% input parameters:
%   x       - vector of parameters [Nn,1]
%   xp      - vector of implicit parameters [Nq,1]
%   f       - frequency vector [Nm,1]
%
% output parameters:
%   R       - response vector only (no frequency), S11 in dB in this case

%% test inputs
% f = 0:0.1e9:20e9;
% x = [5.003 0.002 5.002 0.0075 5.001 0.0075 5 0.03 0.7 1];
% xp = [0.001e-9 5];

%% Scale inputs to be the same as CST
f = f.*1e9;
x([2,4,6]) = x([2,4,6])./1000;

%% Decode inputs
Z0 = 50;
L0 = 2e-3;
Z1 = x(1);
L1 = x(2);
Z2 = x(3);
L2 = x(4);
Zsys = x(5);
LambdaC = x(6);
RatioR = x(7);
% RatioZ = x(8);

IndPin = xp(1);
Zcomp = xp(2);

%% Fixed parameters
N = 10;
PeripheralCoaxIR = 0.62e-3;
CentralCoaxOR = 3.5e-3;
theta2deg = 90;
Lper = 4e-3;
Zper = 5;

%% Get circuit elements
lam = 299792458./f;
bet = 2.*pi./lam;
% Central input coax
Tcoax0 = TlineABCD(Z0,bet,L0);
% Step 01
RI0 = CentralCoaxOR/exp(Z0/60);
RI1 = CentralCoaxOR/exp(Z1/60);
Tstep01 = coaxInStepABCD(2*RI0,2*RI1,2*CentralCoaxOR,f);
% Input coax1
Tcoax1 = TlineABCD(Z1,bet,L1);
% Step 12
RI2 = CentralCoaxOR/exp(Z2/60);
Tstep12 = coaxInStepABCD(2*RI1,2*RI2,2*CentralCoaxOR,f);
% Input coax2
Tcoax2 = TlineABCD(Z2,bet,L2);
% Step 2cone
RIcone = CentralCoaxOR/exp(Zsys/60);
Tstep2cone = coaxInStepABCD(2*RI2,2*RIcone,2*CentralCoaxOR,f);
% Transition
Lt = calcTransitionLength(theta2deg, CentralCoaxOR, Zsys);
Tt = TlineABCD(Zsys,bet,Lt);
% Conical line
% Zper = Zsys*RatioZ;
PeripheralCoaxOR = PeripheralCoaxIR*exp(Zper/60);
Lcone = calcConicalLength(theta2deg, Zsys, CentralCoaxOR, PeripheralCoaxOR, RatioR, N, LambdaC);
Tcone = TlineABCD(Zsys,bet,Lcone);
% Compensating lines
Lcomp = PeripheralCoaxOR*pi/2;
Tcomp = TlineABCD(Zcomp,bet,Lcomp);
% Back short
Lbs = LambdaC/4 - Lcomp;
Tbs = stubABCD(Zsys,bet,Lbs,'p','s'); 
% Series inductor
Tind = LserABCD(IndPin,f);
% Peripheral coax
Tper = TlineABCD(Zper,bet,Lper);

%% Build circuit
[T1,T2t,T2,T3,T] = deal(zeros(2,2,length(f)));
for ff = 1:length(f)
    T1(:,:,ff) = Tcoax0(:,:,ff)*Tstep01(:,:,ff)*Tcoax1(:,:,ff)*Tstep12(:,:,ff)*Tcoax2(:,:,ff)*Tstep2cone(:,:,ff)*...
        Tt(:,:,ff)*Tcone(:,:,ff)*Tcomp(:,:,ff);
    T2t(:,:,ff) = Tcomp(:,:,ff)*Tbs(:,:,ff);
    ZT2 = ABCD2Z(T2t(:,:,ff));
    Zshunt = ZT2(1,1);
    T2(:,:,ff) = ZparABCD(Zshunt);
    T3(:,:,ff) = Tind(:,:,ff)*Tper(:,:,ff);
    T(:,:,ff) = T1(:,:,ff)*T2(:,:,ff)*T3(:,:,ff);
end
Smat = ABCD2S(T,Z0,Zper);
R = dB20(squeeze(Smat(1,1,:)));
% R = squeeze(Smat(1,1,:));

end

function Len = calcTransitionLength(theta2, R2, Z)
    theta2 = theta2*pi/180;
    theta1 = 2*atan(tan(theta2/2)/(exp(2*pi*Z/(120*pi))));
    R1 = R2/exp(Z/60);
    r1 = 3.5*(R2 - R1);
    r2 = (-r1*cos(pi/2 - theta1) - R1*sin(pi/2 - theta1))/(sin(pi/2 - theta1) - 1);
    r0 = (r1 + r2)/2;
    Len = r0*(pi/2 - (pi/2 - theta1)/2);
end

function L = calcConicalLength(theta2, Z, R2, or, RatioR, N, LambdaC)
    Rp = calcPlacementRadius(RatioR, N, LambdaC);
    theta2 = theta2*pi/180;
    theta1 = 2*atan(tan(theta2/2)/(exp(2*pi*Z/(120*pi))));
    R1 = R2/exp(Z/60);
    r1 = 3.5*(R2 - R1);
    r2 = (-r1*cos(pi/2 - theta1) - R1*sin(pi/2 - theta1))/(sin(pi/2 - theta1) - 1);
    rho1 = R2 + r1;
    rho2 = R1 + r2 - r2*sin(pi/2 - theta1);
    lcut = (rho1 + rho2)/2;
    lp = Rp/(cos((pi/2 - theta1)/2)) - lcut;
    L = lp - or*pi/2;
end

function Rp = calcPlacementRadius(RatioR, N, LambdaC)
    Rp = RatioR*(LambdaC/4)*N/pi;
end


