% Script to run the SMmain function
% Use as an example for other cases and to develop the function code

close all
clear all
format compact

%% Set up the models and SM parameters

% Frequencies in GHz
fmin = 1e-9;
fmax = 20;  
Nf = 101;

% Initial input parameters 
xinit = [30 7.5 12 7.5 5 29 0.65]'; 
% xpinit = [0.1e-9 5]';   % Initial implicit parameters 
xpinit = [33e-12 5.05]';   % Initial implicit parameters 
% xpinit = [50e-12 5.00]';   % Initial implicit parameters 


% Set up fine model
Mf.path = 'c:\Work\2015\CST\ICEAA\CombinerSM\';
Mf.name = 'InductanceExtract';
Mf.solver = 'CST';
Mf.params = {'aaa_Z1';'aaa_L1';'aaa_Z2';'aaa_L2';'aaa_Zsys';'aaa_lambda_c';'aaa_ratioR'};
Mf.ximin = [26 5 6 5 4 20 0.6]';
Mf.ximax = [45 10 15 10 7.5 32 0.7]';
% Mf.ximin = [10 5 5.5 5 4.8 25 0.6 0.9]';
% Mf.ximax = [19 10 9 10 5.2 35 0.9 1.1]';
Mf.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);
    
% Set up coarse model
Mc.path = [pwd,'\'];
Mc.name = @conicalCoarseModelDdV;  % Must pass a function handle if MATLAB is the simulator...
Mc.solver = 'MATLAB';
Mc.params = {'x';'xp';'f'}; % Must be this order (and names) for MATLAB sims.  Can omit xp or f though...
Mc.ximin = Mf.ximin;
Mc.ximax = Mf.ximax;
Mc.xpmin = [20e-12 5]';
Mc.xpmax = [200e-12 5.5]';
Mc.freq = reshape(linspace(fmin,fmax,Nf),Nf,1);

% Set up the SM
% The initial SM structure
Sinit.xp = xpinit;

% All the standard SM options - not all shown here... (buildSurr.m for details)
SMopts.getA = 0;
SMopts.getB = 0;
SMopts.getc = 1;
SMopts.getG = 1;
SMopts.getxp = 1;
SMopts.getF = 0;
SMopts.getE = 0;
SMopts.getd = 0;

SMopts.ximin = Mc.ximin;
SMopts.ximax = Mc.ximax;
SMopts.xpmin = Mc.xpmin;
SMopts.xpmax = Mc.xpmax;
SMopts.optsFminS = optimset('display','iter');
SMopts.optsPBIL.display =  'iter';
SMopts.optsPBIL.Nfeval = 5000;
SMopts.wk = 1.15;

% Set up the optimization
OPTopts.Ni = 2;
OPTopts.Rtype = {'S11dB'};
OPTopts.globOpt = 1;
OPTopts.globOptSM = 2;
OPTopts.goalType = {'lt'};
OPTopts.goalResType = {'S11dB'};
OPTopts.goalVal = {-20};
OPTopts.goalWeight = {1};
OPTopts.goalStart = {9};
OPTopts.goalStop = {11};
OPTopts.errNorm = {'L1'};
OPTopts.optsPBIL.display =  'iter'; 
OPTopts.optsPBIL.Nfeval = 5000;
OPTopts.optsPBIL.Nbest = 10; % DOM
OPTopts.M_PBIL = 6;
OPTopts.optsFminS = optimset('display','iter');

%% Run the main loop
[Ri,Si,Pi] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts);

%% Plot some results

