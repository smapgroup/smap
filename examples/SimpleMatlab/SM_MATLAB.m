% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies 
fmin = pi;
fmax = 2*pi;  
Nf = 101;
% f = linspace(0,pi,101);

% Initial input parameters 
xinit = [2,pi/2]';  % Initial input parameters (only ls in this case)
xpinit = [2.1]';   % Initial implicit parameters (only eps_r in this case)


% Set up fine model
Mf.path = [pwd,'\'];
Mf.name = @fineTest;
Mf.solver = 'MATLAB';
Mf.params = {'x','f'};
Mf.ximin = [1,0.5*pi/2]';
Mf.ximax = [3,2*pi/2]';
Mf.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);
    
% Set up coarse model (MATLAB)
Mc.path = [pwd,'\'];
Mc.name = @coarseTest;  % Must pass a function handle if MATLAB is the simulator...
Mc.solver = 'MATLAB';
Mc.params = {'x','f'}; % Must be this order (and names) for MATLAB sims.  Can omit xp or f though...
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = [2]';
Mc.xpmax = [2.2]';
Mc.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);

% Set up the SM
% The initial SM structure
Sinit = [];
% Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 1;
SMopts.getB = 0;
SMopts.getc = 1;
SMopts.getG = 0;
SMopts.getxp = 0;
SMopts.getF = 0;
SMopts.getE = 0;
SMopts.getd = 0;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
% SMopts.xpmin = Mc.xpmin;
% SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','iter');
SMopts.optsPBIL.display =  'iter';
SMopts.optsPBIL.Nfeval = 5000;
% SMopts.wk = 1.15;

% Set up the optimization
OPTopts.Ni = 1;
OPTopts.Rtype = {'Gen'};
OPTopts.globOpt = 0;
OPTopts.globOptSM = 1;
OPTopts.goalType = {'minimax'};
OPTopts.goalResType = {'Gen'};
% OPTopts.goalVal = {-20};
OPTopts.goalWeight = {1};
% OPTopts.goalStart = {1.4e9};
% OPTopts.goalStop = {1.6e9};
OPTopts.errNorm = {'L1'};
OPTopts.optsPBIL.display =  'iter'; 
OPTopts.optsPBIL.Nfeval = 5000;
OPTopts.optsPBIL.Nbest = 10; % DOM
OPTopts.M_PBIL = 4;
OPTopts.optsFminS = optimset('display','iter');
% OPTopts.optsFminS = optimset('MaxFunEvals',10,'display','iter');

%% Run the main loop
[Ri,Si,Pi,Ci,Li] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

%% Plot some results

