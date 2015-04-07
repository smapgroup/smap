function R = MSstubCoarse(x,xp,f)
% input parameters:
%   x       - vector of parameters [Nn,1]
%   xp      - vector of implicit parameters [Nq,1]
%   f       - frequency vector [Nm,1]
%
% output parameters:
%   R       - response vector only (no frequency), S11 in dB in this case

%% test inputs
% f = linspace(1e9,2e9,41);
% x = [80e-9];
% xp = [2.1];

%% Scale inputs to be the same as FEKO
f = f.*1e0;

%% Decode inputs
ls = x(1);

eps_r = xp(1);

%% Fixed parameters
l = 80e-3;
w = 5e-3;
h = 1.5e-3;
t = 0;

%% Get circuit elements
% Get ustrip impedance and propagation constant
[Zf,bet] = ustripLine(eps_r,h,w,t,f);

% FeedLine
Tfeed = TlineABCD(Zf,bet,l/2);
% Stub
Tstub = TlineABCD(Zf,bet,ls);

%% Build circuit
[Ts,T] = deal(zeros(2,2,length(f)));
for ff = 1:length(f)
    ZTs = ABCD2Z(Tstub(:,:,ff));
    Zshunt = ZTs(1,1);
    Ts(:,:,ff) = ZparABCD(Zshunt);
    T(:,:,ff) = Tfeed(:,:,ff)*Ts(:,:,ff)*Tfeed(:,:,ff);
end
Smat = ABCD2S(T,50,50);
R = dB20(squeeze(Smat(1,1,:)));
% R = squeeze(Smat(1,1,:));

end


