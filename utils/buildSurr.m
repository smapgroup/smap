function Si = buildSurr(xi,Rfi,S,opts)

% function Si = buildSurr(xi,Rfi,S,opts)
% function to build a surrogate model structure S, using the input fine 
% model response pairs {xi,Rfi}.  More than one input response pair can be
% specified, as is typical in space mapping modelling instead of
% optimization applications.  User specified options - mainly flags to
% specify SM types - are provided in opts.
%
% Based on the 2006 MTT paper by Koziel et al on a SM Framework
%
% In general the extraction performs best when all input variables,
% including implicit variables, are of the degree 1.  This is because the
% internal optimizer attemps to optimize additive as well as multiplicative
% variables - the latter typically normalized to 1 - at the same time.
%
% Inputs:
% xi: Response positions {1,Ni}[Nn,1] (Nn input parameters, Ni previous fine evaluations)
% Rfi: Fine model responses {1,Ni}[Nm,1] (Nm output responses, Ni previous fine evaluations)
%   xi and Rfi may be cell arrays of the same length indicating multiple fine
%   model evaluations.  The surrogate will then be a best fit over the entire
%   extended parameter space.
% S:  Initial surrogate model containing at least
%  S.coarse -  Function handle pointing to a predefined function that evaluates
%          the coarse model operating as Rc = coarse(xc,xp,f) thus returning Rc 
%          (coarse model response) for each element in xc. Optional xp are
%          the implicit (pre-assigned) parameters.
%  S.xp - Optional pre-assigned variables (depending on the coarse model [Nq,1]
%  S.f - Optional frequency vector where the responses are calculated [Nm,1] (must be included for FSM)
% opts: User options structure - mostly flags to specifiy type of SM (defaults)
%   getA - flag to get multiplicative OSM (set to 2 to get M factors), typically 1 to force single factor for all outputs which is much faster and often works well (0)
%   getB - fleg to get multiplicative input SM (0)
%   getc - flag to get additive input SM (1)
%   getG - flag to get multiplicative ISM (0)
%   getxp- fleg to optimize additive ISM (0)
%   getd - flag to get additive OSM (0)
%   getE - flag to get first order OSM (0)
%   getF - flag to get frequency mapping (0)
%   wk - weights to determine the error function in the model fitting when
%        more than one fine point is included.  Can be the same length 
%        as the number of cells in xi and Rfi.  Default all (1). See eq
%        (11) in the reference. If length(xi) > 1 and length(wk) == 1 then
%        wk will be generated internally as an exponential function,
%        growing by the factor wk ie: wk.^[1:length(xi)]
%   vk - weights to determine the error function in the Jacobian fitting when
%        more than one fine point is included.  Can be the same length 
%        as the number of cells in xi and Rfi.  Default all (0). See eq
%        (11) in the reference. If length(xi) > 1 and length(vk) == 1 then
%        wk will be generated internally as an exponential function,
%        growing by the factor vk ie: vk.^[1:length(xi)]
%   ximin - vector of minimum values for xi to constrain the search space
%   ximax - vector of maximum values for xi to constrain the search space
%   xpmin - vector of minimum values for xp to constrain the search space
%   xpmax - vector of maximum values for xp to constrain the search space
%   AlimMin - minimum value for A to constrain the search space
%   AlimMax - maximum value for A to constrain the search space
%   Fmin - minimum values for F to constrain the search space [2,1]
%   Fmax - maximum values for F to constrain the search space [2,1]
%   optsFminS - options for the fminsearch local optimizer in the PE
%   globOpt - flag to include global search (PBILreal) in the PE 
%   M_PBIL - number of bits per parameter in the PBIL global search (8)
%   optsPBIL - options for the PBIL global optimizer in the PE
%   errNorm - Type of error norm to use for parameter extraction ('L1','L2')
%             Default is the L1 norm (Least absolute deviations)
%   errW - Vector of weights (typically binary but can be any real number),
%          of length Nm, to calculate the extraction error.  Can be used to 
%          mask out regions in the response domain of less importance. 
%          Default ones. 
%
% Returns:
% S -  The surrogate model structure containing:
% coarse:  Function handle pointing to a predefined function that evaluates
%          the coarse model operating as Rc = coarse(xc,xp) thus returning Rc 
%          (coarse model response) for each element in xc. Optional xp are
%          the implicit (pre-assigned) parameters.
% A:       Multiplicative OSM factor diag[Nm,Nm]
% B:       Multiplicative input SM factor [Nn,Nn]
% c:       Additative input SM term [Nn,1]
% G:       Multiplicative ISM factor [Nq,Nm]
% xp:      Pre-assigned parameters [Nq,1]
% d:       Additive zeroth order OSM term [Nm,1]
% E:       Additive first order OSM term [Nm,Nn]
% xi:      Position of last update point of the model [Nn,1]
% F:       Frequency space mapping parameters [2,1]
% freq:    Coarse model frequency range [Nm,1]
%
% Date created: 2014-11-09
% Dirk de Villiers
% Last Modified: 2015-03-12
% Updates:
% 2014-11-09: Write function shell and basic functionality
% 2014-11-10: Fix constraints of input SM (case 1), cases 2 - 5
% 2014-11-11: Impliment add OSM (d), Fix B to diagonal matrix
% 2014-11-13: Add some more 2 variable cases. Include global search flag
%             and constraint flags. Add global search to case 5.
% 2014-11-14: Add error norm and error weight options
% 2014-11-22: Move some general code out of switch
% 2014-11-23: Include optimizer options as inputs
% 2015-02-23: Fix implicit constraints
% 2015-03-03: Fix simple special case of only getA for one fine model
%             Include c and B cases, and simplify Bc case (no constraints)   
% 2015-03-04: Include AB, Ac, ABc, G, xp, cxp cases
%             Simplify Gxp case (no constraints)
%             Include AlimMin/Max as opts parameter 
%             Include Gv_init variable for consistency with Bv_init
% 2015-03-05: Include globalOpt in BcGxp case
%             Include BG, cG, BcG, BGxp, cGxp cases
%             REDESIGN CASES COMPLETELY.  ALL GENERAL.
% 2015-03-09: Start FSM
% 2015-03-10: Continue FSM
% 2015-03-11: Finish FSM
% 2015-03-12: Constrained search added
% 2015-03-13: Single value wk functionality added
% ToDo: Impliment E (first order OSM)
% ToDo: More optimizer options - just use fminsearch now
% ToDo: Jacobian fitting in error functions (vk)

% Preassign some variables
[lenA,lenB,lenc,lenG,lenxp,lenF] = deal(0);
[Nc,Nn,Nq,Nm] = deal(0);
[Av_init,Bv_init,c_init,Gv_init,xp_init,F_init,typeA,typeB,typec,typeG,typexp,typeF] = deal([]);
[Amin,Amax,Bmin,Bmax,cmin,cmax,Gmin,Gmax,xpmin,xpmax,Fmin,Fmax] = deal([]);

% Sort out formats
if iscell(xi) && iscell(Rfi) && (length(xi) == length(Rfi))   % Basic error checking - if any issue here only first point will be used
    Nc = length(xi);  % Number of input point cells
else     % Force only first point to be used, and make cell arrays
    if ~iscell(xi), xi = mat2cell(xi); end 
    if ~iscell(Rfi), Rfi = mat2cell(Rfi); end 
    Nc = 1;
end

% Set defaults and do BASIC error checking
getA = 0;
getB = 0;
getc = 1;
getG = 0;
getxp = 0;
getd = 0;
getE = 0;
getF = 0;
ximin = 0.8.*xi{Nc}; % Default constraints 
ximax = xi{Nc}./0.8;
if isfield(S,'xp')
    xpmin(:,1) = 0.8.*S.xp; % Default constraints 
    xpmax(:,1) = S.xp./0.8;
    Nq = length(S.xp);
end
AlimMin = 0.5;  
AlimMax = 2;
if isfield(S,'f')
    FminDef = [0.5;-mean(S.f)];    % Default constraints
    FmaxDef = [2;mean(S.f)];
end
optsFminS = optimset('display','none');
globOpt = 0;
M_PBIL = 8;
optsPBIL = [];
errNorm = 'L1';
errW = 1;
if isfield(opts,'wk') 
    wk = opts.wk;
    if length(wk) == 1
        wk = wk.^[1:Nc];
    end
else
    wk = ones(1,Nc);
end
if isfield(opts,'vk') 
    vk = opts.vk;
    if length(vk) == 1
        vk = vk.^[1:Nc];
    end
else
    vk = ones(1,Nc);
end

if isfield(opts,'getA'), getA = opts.getA; end
if isfield(opts,'getB'), getB = opts.getB; end
if isfield(opts,'getc'), getc = opts.getc; end
if isfield(opts,'getG'), getG = opts.getG; end
if isfield(opts,'getxp'), getxp = opts.getxp; end
if isfield(opts,'getd'), getd = opts.getd; end
if isfield(opts,'getE'), getE = opts.getE; end
if isfield(opts,'getF'), getF = opts.getF; end
if isfield(opts,'ximin'), ximin = opts.ximin; end
if isfield(opts,'ximax'), ximax = opts.ximax; end
if isfield(opts,'xpmin'), xpmin = opts.xpmin; end
if isfield(opts,'xpmax'), xpmax = opts.xpmax; end
if isfield(opts,'AlimMin'), AlimMin = opts.AlimMin; end
if isfield(opts,'AlimMax'), AlimMax = opts.AlimMax; end
if isfield(opts,'Fmin'), FminDef = opts.Fmin; end
if isfield(opts,'Fmax'), FmaxDef = opts.Fmax; end
if isfield(opts,'optsFminS'), optsFminS = opts.optsFminS; end
if isfield(opts,'globOpt'), globOpt = opts.globOpt; end
if isfield(opts,'M_PBIL'), M_PBIL = opts.M_PBIL; end
if isfield(opts,'optsPBIL'), optsPBIL = opts.optsPBIL; end
if isfield(opts,'errNorm'), errNorm = opts.errNorm; end
if isfield(opts,'errW'), errW = opts.errW; end

% Limit the parameter space
S.ximin = ximin;
S.ximax = ximax;
if isfield(S,'xp')
    S.xpmin = xpmin;
    S.xpmax = xpmax;
end

% Get vector sizes
Nn = length(xi{Nc});
Si = S;



% Initialise the optimization variable vector
% The limits named *min/*max indicate the global search range.  Actual
% values may fall outside these ranges after the local search alignment...
if getA
    typeA = 'A';
    if isfield(S,'A')
        A_init = S.A;   % Try to always initialize - avoid running coarse model for no real reason...
        Nm = length(diag(A_init));
        Av_init = diag(A_init);
        Amin = Av_init.*AlimMin;
        Amax = Av_init.*AlimMax;
        % else case only handled if needed
    else
        % Have to run the coarse model to find the response size
        Rc = evalSurr(xi{Nc},S);
        [Nm,Np] = size(Rc);
        A_init = eye(Nm);
        Av_init = diag(A_init);
        Amin = Av_init.*AlimMin;
        Amax = Av_init.*AlimMax;
        lenA = Nm;
    end
    if getA == 1 % Use a single A for all responses
        optsParE.A = mean(Av_init);
        optsParE.Nm = Nm;
        Av_init = mean(Av_init);
        Amin = Av_init.*AlimMin;
        Amax = Av_init.*AlimMax;
        lenA = 1;
    end
end
if getB
    typeB = 'B';
    lenB = Nn;
    if isfield(S,'B')
        B_init = S.B;
    else
        B_init = eye(Nn);
    end
    Bv_init = diag(B_init);
    Bmin = ximin./xi{Nc};    % Assume c = 0
    Bmax = ximax./xi{Nc};
end
if getc
    typec = 'c';
    lenc = Nn; 
    if isfield(S,'c')
        c_init(:,1) = S.c;
    else
        c_init = zeros(Nn,1);
    end
    cmin = ximin - xi{Nc};   % Assume B = 1
    cmax = ximax - xi{Nc};
end
% Have to provide an xp input for implicit space mapping...
if getxp && isfield(S,'xp')
    typexp = 'xp';
    lenxp = Nq; 
    xp_init(:,1) = S.xp;
    Nq = length(xp_init);
end
if getG && isfield(S,'xp')
    typeG = 'G';
    lenG = Nq*Nn;
    if isfield(S,'G')
        G_init = S.G;
    else
        G_init = zeros(Nq,Nn);  % In case there is no 'G'
    end
    Gv_init = [reshape(G_init,Nq*Nn,1)];
    % Estimate edge values here per row - system underdetermined...
    % Assume xp = 0, and hope everything is well scaled!  Otherwise the
    % local search will have to sort out the issues.
    xsum = sum(xi{Nc});
    [GminM,GmaxM] = deal(ones(Nq,Nn));
    for qq = 1:Nq
        GminM(qq,:) = ones(1,Nn)*xpmin(qq)/xsum;
        GmaxM(qq,:) = ones(1,Nn)*xpmax(qq)/xsum;
    end
    Gmin = reshape(GminM,Nq*Nn,1);
    Gmax = reshape(GmaxM,Nq*Nn,1);
end
% Force xpmin/xpmax to empty if xp not requested
if ~getxp, [xpmin,xpmax] = deal([]); end
if getF && isfield(S,'f')
    typeF = 'F';
    lenF = 2;
    if isfield(S,'F')
        F_init = S.F;
    else
        F_init = [1;0];
    end
    Fmin = FminDef; % Take defaults directly
    Fmax = FmaxDef;
end
    
initVect = [Av_init;Bv_init;c_init;Gv_init;xp_init;F_init];
minVect = [Amin;Bmin;cmin;Gmin;xpmin;Fmin];
maxVect = [Amax;Bmax;cmax;Gmax;xpmax;Fmax];
inputType = [typeA,typeB,typec,typeG,typexp,typeF];
Si.inputType = inputType;

% Set up the general parameter extraction options
optsParE.Nn = Nn;
optsParE.Nq = Nq;
optsParE.lenA = lenA;
optsParE.lenB = lenB;
optsParE.lenc = lenc;
optsParE.lenG = lenG;
optsParE.lenxp = lenxp;
optsParE.lenF = lenF;
optsParE.errNorm = errNorm;
optsParE.errW = errW;

% Set up and run the optimizations (parameter extractions)
% Use provided SM parameters as initial values if available, and enhance
% with global search if required
if strcmp(inputType,'F') || strcmp(inputType,'AF') && Nc == 1 % Special cases where the coarse model is not re-evaluated (for each optimization iteration).  Use interpolation/extrapolation instead...
    if ~exist('Rc','var')    % Check if coarse model has been calculated
        Rc = evalSurr(xi{Nc},S);
    end
    if globOpt
        F_init = PBILreal(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),Fmin,Fmax,M_PBIL,optsPBIL);
    end
    Fvect = fminsearch(@(Fvect) erriF(Fvect,Rfi,Rc,S.f,optsParE),F_init,optsFminS);
    fs = Fvect(1).*S.f + Fvect(2);
    Rs = interp1(S.f,Rc,fs,'spline');
    if strcmp(inputType,'AF')       % Need to also get the A factor
        A = Rfi{1}./reshape(Rs,Nm,1);
        if getA == 1, A = mean(A); end
        optVect = [A;reshape(Fvect,lenF,1)];
    else
        optVect = reshape(Fvect,lenF,1);
    end
elseif strcmp(inputType,'A') && Nc == 1  % Special case without optimization
    if ~exist('Rc','var')    % Check if coarse model has been calculated
        Rc = evalSurr(xi{Nc},S);
    end
    A = Rfi{1}./Rc;
    if getA == 1, A = mean(A); end
    optVect = A;
else
    if globOpt
        % Start with global search to get initial value
        [initVect] = PBILreal(@(optVect) erri(optVect,xi,Rfi,S,wk,vk,optsParE),minVect,maxVect,M_PBIL,optsPBIL);
    end
    optVect = fminsearch(@(optVect) erri(optVect,xi,Rfi,S,wk,vk,optsParE),initVect,optsFminS);
end
% Rebuild the individual parameters from the vector
A = diag(optVect(1:lenA));
B = diag(optVect(1+lenA:lenA+lenB));
c = reshape(optVect(1+lenA+lenB:lenA+lenB+lenc),lenc,1);
G = reshape(optVect(1+lenA+lenB+lenc:lenA+lenB+lenc+lenG),min(lenG,Nq),Nn); % Must be empty matrix if lenG == 0
xp = reshape(optVect(1+lenA+lenB+lenc+lenG:lenA+lenB+lenc+lenG+lenxp),lenxp,1);
F = reshape(optVect(1+lenA+lenB+lenc+lenG+lenxp:lenA+lenB+lenc+lenG+lenxp+lenF),lenF,1); % Must be empty matrix if lenG == 0

if lenA > 0, Si.A = A; end
if lenB > 0, Si.B = B; end
if lenc > 0, Si.c = c; end
if lenG > 0, Si.G = G; end
if lenxp > 0, Si.xp = xp; end
if lenF > 0, Si.F = F; end

% Additive zero order OSM 
if getd, Si.d = Rfi{end} - evalSurr(xi{end},Si); end

% Additive first order OSM - ToDo

end


% Error function for optimization
function e = erri(optVect,xi,Rfi,S,wk,vk,opts)
% Extract individual parameters
A = diag(optVect(1:opts.lenA));
B = diag(optVect(1+opts.lenA:opts.lenA+opts.lenB));
c = reshape(optVect(1+opts.lenA+opts.lenB:opts.lenA+opts.lenB+opts.lenc),opts.lenc,1);
G = reshape(optVect(1+opts.lenA+opts.lenB+opts.lenc:opts.lenA+opts.lenB+opts.lenc+opts.lenG),min(opts.lenG,opts.Nq),opts.Nn);    % Must be empty matrix if lenG == 0...
xp = reshape(optVect(1+opts.lenA+opts.lenB+opts.lenc+opts.lenG:opts.lenA+opts.lenB+opts.lenc+opts.lenG+opts.lenxp),opts.lenxp,1);
F = reshape(optVect(1+opts.lenA+opts.lenB+opts.lenc+opts.lenG+opts.lenxp:opts.lenA+opts.lenB+opts.lenc+opts.lenG+opts.lenxp+opts.lenF),opts.lenF,1);
% Update surrogate model structure
if opts.lenA > 0, S.A = A; end
if opts.lenB > 0, S.B = B; end
if opts.lenc > 0, S.c = c; end
if opts.lenG > 0, S.G = G; end
if opts.lenxp > 0, S.xp = xp; end
if opts.lenF > 0, S.F = F; end
% Calculate the error function value
Nc = length(wk);
ec = zeros(1,Nc);
if length(opts.errW) == 1
    errW = Rfi{1}./Rfi{1};
else
    errW = opts.errW;
end
for cc = 1:Nc
    Rs = evalSurr(xi{cc},S);
    diffR = errW.*(Rfi{cc} - Rs);
    if isequal(opts.errNorm,'L1')
        ev = sum(abs(diffR),2);  % Error vector [Nm,1]
    elseif isequal(opts.errNorm,'L2')
        ev = sum(abs(diffR).^2,2);  % Error vector [Nm,1]
    end
    ec(cc) = wk(cc).*sum(ev);    % This can be updated if more specific error functions are needed as function of response
end
e = sum(ec);
end


% Special case error function where only F is optimized and the coarse
% model is not re-evaluated - interpolation/extrapolation is used on the
% provided coarse model response...
function e = erriF(Fvect,Rfi,Rc,f,opts)

fs = Fvect(1).*f + Fvect(2);
Rs = interp1(f,Rc,fs,'spline');

diffR = Rfi{1} - reshape(Rs,length(f),1);
if isequal(opts.errNorm,'L1')
    e = sum(abs(diffR));  % Error vector [Nm,1]
elseif isequal(opts.errNorm,'L2')
    e = sum(abs(diffR).^2);  % Error vector [Nm,1]
end
end

