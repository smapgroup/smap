function [Ri,Si,Pi,Ci, Li] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts)

% Space Mapping main loop

% Inputs:
% xinit:    Initial values of input parameters [Nn,1]
% Sinit:    Initial SM structure (see buildSurr.m/evalSurr.m)
%           Important to include initial implicit parameters
% SMopts:   SM options (see buildSurr.m/evalSurr.m)
% Mf:       Fine model structure
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% Mc:       Coarse model structure (can be cell array of structures if more than one type has to be calculated to get all the fine model responses)
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% OPTopts:  Optimization options structure for external loop
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types (so far):
%               'S11dB' 
%               'S11complex'
%   Ni:         Maximum number of iterations
%   globOpt:    Flag to run PBIL (1 for only first iteration, 2 for all iterations) (default 0)
%   M_PBIL:     Vector of bits for the global search variables (see PBILreal.m)
%   globOptSM:  Flag to run PBIL during the PE process (1 for only first iteration, 2 for all iterations) (default 0)
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11dB'
%   goalType:   Cell array of goal types {1,Ng}
%               Valid types:
%                   'lt' (Less than)
%                   'gt' (Greater than)
%                   'eq' - todo
%                   'minimax' - todo
%   goalVal:    Cell array of goal values {1,Ng} - same order as goalType
%   goalWeight: Vector of goal weights [1,Ng]
%   goalStart:  Cell array of start of valid goal domain {1,Ng}
%   goalStop:   Cell array of stop of valid goal domain {1,Ng}
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)
%   optsFminS:  options for the fminsearch local optimizer
%   optsPBIL:   options for the PBIL global optimizer
%   plotIter:   Flag to plot the responses after each iteration


% Returns:
% Ri:   Structure containing the responses at each iteration
% Si:   Structure containing the surrogates at each iteration
% Pi:   Structure containing the parameters at each iteration
% Ci:   Structure containing the costs at each iteration
% Li:   Structure containing the limiting information at each iteration

% Date created: 2015-03-06
% Dirk de Villiers and Ryno Beyers
% Last Modified: 2015-03-27
% Updates:
% 2015-03-06: Write function shell and basic functionality
% 2015-03-09: Continue with shell and basic functionality
% 2015-03-10: Continue with shell and basic functionality
% 2015-03-11: Continue with shell and basic functionality
% 2015-03-12: Continue with shell and basic functionality
% 2015-03-13: Continue with shell and basic functionality
% 2015-03-14: Add MATLAB case to fine model evaluation function
% 2015-03-27: Add costeq cost function
%             Add plotIter functionality


% Set defaults
Ni = 10;    % Maximum unmber of iterations
globOpt = 0;
M_PBIL = 8;
globOptSM = 0;
optsFminS = optimset('display','none');
optsPBIL = [];
plotIter = 1;

if isfield(OPTopts,'Ni'), Ni = OPTopts.Ni; end
if isfield(OPTopts,'globOpt'), globOpt = OPTopts.globOpt; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'globOptSM'), globOptSM = OPTopts.globOptSM; end
if isfield(OPTopts,'optsFminS'), optsFminS = OPTopts.optsFminS; end
if isfield(OPTopts,'M_PBIL'), M_PBIL = OPTopts.M_PBIL; end
if isfield(OPTopts,'optsPBIL'), optsPBIL = OPTopts.optsPBIL; end
if isfield(OPTopts,'plotIter'), plotIter = OPTopts.plotIter; end

% Set up models - bookkeeping
Nq = 0;
Nn = length(xinit);
xi{1} = reshape(xinit,Nn,1);
if isfield(Sinit,'xp')
    Nq = length(Sinit.xp);
    xpi{1} = reshape(Sinit.xp,Nq,1);
end
Nr = length(OPTopts.Rtype); % Number of responses requested
Mc.Rtype = OPTopts.Rtype;
Mf.Rtype = OPTopts.Rtype;
Sinit.M = Mc;
Sinit.coarse = @coarseMod;
if isfield(Mc,'freq')
    Sinit.f = Mc.freq;
    fc = Mc.freq;
else
    fc = [];
end
useAllFine = 0;
if isfield(SMopts,'wk')
    % Assume user wants to use all calculated fine responses to fit the
    % surrogate...
    useAllFine = 1;
end

% Optimize coarse model to find initial alignment position
% Limit the parameter space
if isfield(Mc,'ximin'), Sinit.ximin = Mc.ximin; end
if isfield(Mc,'ximax'), Sinit.ximax = Mc.ximax; end
if isfield(Mc,'xpmin'), Sinit.xpmin = Mc.xpmin; end
if isfield(Mc,'xpmax'), Sinit.xpmax = Mc.xpmax; end

if globOpt
    [xinit,costSi,exitFlag,output] = PBILreal(@(xi) costSurr(xi,Sinit,OPTopts),Mc.ximin,Mc.ximax,M_PBIL,optsPBIL);
    xinit = reshape(xinit,Nn,1);
end
[xi{1}, costSi] = fminsearch(@(xi) costSurr(xi,Sinit,OPTopts),xinit,optsFminS);

% Enter the main loop
% [Rci,Rfi,Rsi,Rsai,Si] = deal(cell(1,Ni));
% Rsai is the aligned surrogate, and Rsi the optimized surrogate at each
% iteration
costC = costSi; % Only calculate the coarse model cost once, at iteration 0
specF = 0;  % Flag to test if the fine model reached spec
ii = 0;
while ii <= Ni && ~specF
    ii = ii+1;
    Rci{ii} = coarseMod(Mc,xi{ii},Sinit.xp,fc);
    Rfi{ii} = fineMod(Mf,xi{ii});
    if ii == 1
        for rr = 1:Nr
            if globOptSM > 0, SMopts.globOpt = 1; end
            Si{ii}{rr} = buildSurr(xi{ii},Rfi{ii}{rr}.r,Sinit,SMopts);
            Rsi{ii}{rr}.r = evalSurr(xi{ii},Si{ii}{rr});
            Rsi{ii}{rr}.t = Rci{ii}{rr}.t;
            if isfield(Rci{ii}{rr},'f'), Rsi{ii}{rr}.f = Rci{ii}{rr}.f; end
            Rsai{ii}{rr} = Rsi{ii}{rr};
        end
    else
        for rr = 1:Nr
            Rsi{ii}{rr}.r = evalSurr(xi{ii},Si{ii-1}{rr});
            Rsi{ii}{rr}.t = Rci{ii}{rr}.t;
            if isfield(Rci{ii}{rr},'f'), Rsi{ii}{rr}.f = Rci{ii}{rr}.f; end
            if globOptSM < 2, SMopts.globOpt = 0; end
            if ~useAllFine
                Si{ii}{rr} = buildSurr(xi{ii},Rfi{ii}{rr}.r,Si{ii-1}{rr},SMopts);
            else
                for iii = 1:ii
                    r{iii} = Rfi{iii}{rr}.r;
                end
                Si{ii}{rr} = buildSurr(xi,r,Si{ii-1}{rr},SMopts);
            end
            % Also get the currently aligned surrogate for comparison
            Rsai{ii}{rr}.r = evalSurr(xi{ii},Si{ii}{rr});
            Rsai{ii}{rr}.t = Rci{ii}{rr}.t;
            if isfield(Rci{ii}{rr},'f'), Rsai{ii}{rr}.f = Rci{ii}{rr}.f; end
        end
    end
    
    % Test fine model response
    costFi = costFunc(Rfi{ii},OPTopts);
    if costFi == 0  % Exit is spec is reached (will typically not work for eq and never for minimax...)
        specF = 1;
    else
        specF = 0;
        % Do optimization
        if (ii == 1 && globOpt) || globOpt == 2
            [xinit,costSi,exitFlag,output] = PBILreal(@(xi) costSurr(xi,Si{ii}{:},OPTopts),Mc.ximin,Mc.ximax,M_PBIL,optsPBIL);
            xinit = reshape(xinit,Nn,1);
        end
        [xi{ii+1}, costSi] = fminsearch(@(xi) costSurr(xi,Si{ii}{:},OPTopts),xinit,optsFminS);
        
        % Check if inputs will be limited, and don't actually limit the optimized
        % position here.  Important since the SM might move it back within
        % limits in the next iteration, therefore changing the evaluation
        % position of the surrogate model and effectively changing the
        % surrogate parameters!
        % This can make the optimizer very inefficient.  Check if a parameter
        % is limited and try to remove it from the optimization space
        % externally by changing the models.
        % Return vector with 0 if not, 1 for max and -1 for min...
        [limF{ii},limC{ii}] = deal(zeros(Nn,1));
        % Fine model
        if isfield(Mc,'ximin')
            minI = xi{ii+1} < Mf.ximin;
            limF{ii}(minI) = -1;
        end
        if isfield(Mc,'ximax')
            maxI = xi{ii+1} > Mf.ximax;
            limF{ii}(maxI) = 1;
        end
        % Coarse model
        if isfield(Mc,'ximin')
            minI = xi{ii+1} < Mc.ximin;
            limC{ii}(minI) = -1;
        end
        if isfield(Mc,'ximax')
            maxI = xi{ii+1} > Mc.ximax;
            limC{ii}(maxI) = 1;
        end
        % Surrogate model (bit more complicated!)
        for rr = 1:Nr
            limS{ii}{rr} = zeros(Nn,1);
            % First get some defaults or actual values from the surrogate...
            if ~isfield(Si{ii}{rr},'B')
                B = eye(Nn);
            else
                B = Si{ii}{rr}.B;
            end
            if ~isfield(Si{ii}{rr},'c')
                c = zeros(Nn,1);
            else
                c = Si{ii}{rr}.c;
            end
            xis = B*xi{ii+1} + c;
            if isfield(Mc,'ximin')
                minI = xis < Mc.ximin;
                limS{ii}{rr}(minI) = -1;
            end
            if isfield(Mc,'ximax')
                maxI = xis > Mc.ximax;
                limS{ii}{rr}(maxI) = 1;
            end
        end
    end
    
    
    costS{ii} = costSi;
    costF{ii} = costFi;
    % Make a (crude) log file
    save SMlog ii xi Rci Rfi Rsi Si costS costF limF limC limS
    
    if plotIter
        figure(ii)
        for rr = 1:Nr
            subplot(Nr,1,rr)
            if isfield(Rci{ii}{rr},'f')
                plot(Rci{ii}{rr}.f,Rfi{ii}{rr}.r,'k'), grid on, hold on
                plot(Rci{ii}{rr}.f,Rci{ii}{rr}.r,'r'), grid on, hold on
                plot(Rsi{ii}{rr}.f,Rsi{ii}{rr}.r,'b')
                plot(Rsai{ii}{rr}.f,Rsai{ii}{rr}.r,'g')
                xlabel('Frequency')
                % Plot the specs...
            else
                plot(Rfi{ii}{rr}.r,'k'), grid on, hold on
                plot(Rci{ii}{rr}.r,'r'), grid on, hold on
                plot(Rsi{ii}{rr}.r,'b')
                plot(Rsai{ii}{rr}.r,'g')
                xlabel('Index')
                % Plot the specs...
            end
            ylabel(OPTopts.Rtype{rr})
            title(['Iteration ',num2str(ii)])
            legend('Fine','Coarse','Optimised Surrogate', 'Aligned Surrogate')
        end
        
    end
    
end


% Handle output structures
Ri.Rc = Rci;
Ri.Rf = Rfi;
Ri.Rs = Rsi;    % Surrogate after optimization
Ri.Rsa = Rsai;  % Surrogate before optimization, just after alignment at end of previous iteration

Pi = xi;

Ci.costC = costC;
Ci.costS = costS;
Ci.costF = costF;

Li.limC = limC;
Li.limS = limS;
Li.limF = limF;

end



function Rf = fineMod(M,xi)

% Rf is a cell array of structures containing the response in Rf.r, the type Rf.t, and the
% (optional) domain (typically frequency) in Rf.f.  Same length as M.Rtype
% xi is an array of input paraeters - same order as those specified in M
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'S11dB' - obvious!

% Limit the inputs
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
end
if isfield(M,'ximax')
    maxI = xi > M.ximax;
    xi(maxI) = M.ximax(maxI);
end

Nn = length(xi);
% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end
Nr = length(Rtype);
Rf = cell(1,Nr);

% If solver is CST:
if strcmp(M.solver,'CST')
    % Start CST activeX
    cst = actxserver('CSTSTUDIO.Application');
    % Get handle to the model
    mws = invoke(cst,'OpenFile',[M.path,M.name,'.cst']);
    % Update parameters
    for nn = 1:Nn
        invoke(mws,'StoreParameter',M.params{nn},xi(nn));
    end
    invoke(mws,'Rebuild');
    % Run simulation
    solver = invoke(mws,'Solver');
    invoke(solver,'Start');
    
    % Generate output
    for rr = 1:Nr
        if strcmp(Rtype{rr},'S11dB')
            result = invoke(mws,'Result1D','d1(1)1(1)');    % S11 in dB
            % Get nr of frequency points in the plot
            nRead = invoke(result,'GetN');
            [fin,S11in] = deal(zeros(nRead,1));
            for nn = 1:nRead
                fin(nn) = invoke(result,'GetX',nn-1);        % Typically in GHz
                S11in(nn) = invoke(result,'GetY',nn-1);
            end
            if isfield(M,'freq')
                Nm = length(M.freq);
                Rf{rr}.r = reshape(interp1(fin,S11in,M.freq,'spline'),Nm,1);
                Rf{rr}.f = M.freq;
            else
                Nm = nRead;
                Rf{rr}.r = S11in;
                Rf{rr}.f = fin;
            end
            Rf{rr}.t = Rtype{rr};
            release(result);
        elseif strcmp(Rtype{rr},'S11complex')
            resultA = invoke(mws,'Result1D','a1(1)1(1)');    % amplitude of S11
            resultP = invoke(mws,'Result1D','p1(1)1(1)');    % amplitude of S11
            % Get nr of frequency points in the plots
            nRead = invoke(resultA,'GetN');
            [fin,S11in] = deal(zeros(nRead,1));
            for nn = 1:nRead
                fin(nn) = invoke(resultA,'GetX',nn-1);        % Typically in GHz
                amp = invoke(resultA,'GetY',nn-1);
                phase = rad(invoke(resultP,'GetY',nn-1));
                S11in(nn) = amp.*exp(1i*phase);
            end
            if isfield(M,'freq')
                Nm = length(M.freq);
                Rreal = reshape(interp1(fin,real(S11in),M.freq,'spline'),Nm,1);
                Rimag = reshape(interp1(fin,imag(S11in),M.freq,'spline'),Nm,1);
                Rf{rr}.r = Rreal + 1i*Rimag;
                Rf{rr}.f = M.freq;
            else
                Nm = nRead;
                Rf{rr}.r = S11in;
                Rf{rr}.f = fin;
            end
            Rf{rr}.t = Rtype{rr};
            release(resultA);
            release(resultP);
        end
    end
    invoke(mws,'Save');
    invoke(mws,'Quit');
elseif strcmp(M.solver,'MATLAB')    % If solver is MATLAB
    Ni = length(M.params);  % This is interpreted as the number of inputs to the function
    inType = [];
    for ii = 1:Ni
        inType = [inType,M.params{ii}];   % Initialise
    end
    switch inType
        case 'xf'
            Rfi = M.name(xi,f);
        otherwise
            Rfi = M.name(xi);
    end
    % Distribute the responses
    % For MATLAB case the model should the return the specified responses
    % columnwise...
    for rr = 1:Nr
        Rf{rr}.r = Rfi(:,rr);
        Rf{rr}.t = Rtype{rr};
        if exist('f','var')
            Rf{rr}.f = f;
        end
    end
end
end

function Rc = coarseMod(M,xi,xp,f)

% Rc is a cell array of structures containing the response in Rc.r, the type Rc.t, and the
% (optional) domain (typically frequency) in Rc.f.  Same length as M.Rtype
% xi is an array of input parameters - same order as those specified in M
% xp is an array of implicit parameters - same order as those specified in M
% f is an array of frequency points where to evaluate the model (optional)
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'S11dB' - obvious!

% Limit the inputs
if isfield(M,'ximin')
    minI = xi < M.ximin;
    xi(minI) = M.ximin(minI);
end
if isfield(M,'ximax')
    maxI = xi > M.ximax;
    xi(maxI) = M.ximax(maxI);
end
if isfield(M,'xpmin')
    minIp = xp < M.xpmin;
    xp(minIp) = M.xpmin(minIp);
end
if isfield(M,'xpmax')
    maxIp = xp > M.xpmax;
    xp(maxIp) = M.xpmax(maxIp);
end

Nn = length(xi);
% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end
Nr = length(Rtype);
Rc = cell(1,Nr);

% If solver is MATLAB:
if strcmp(M.solver,'MATLAB')
    Ni = length(M.params);  % This is interpreted as the number of inputs to the function
    inType = [];
    for ii = 1:Ni
        inType = [inType,M.params{ii}];   % Initialise
    end
    switch inType
        case 'xxpf'
            Rci = M.name(xi,xp,f);
        case 'xf'
            Rci = M.name(xi,f);
        case 'xxp'
            Rci = M.name(xi,xp);
        otherwise
            Rci = M.name(xi);
    end
    % Distribute the responses
    % For MATLAB case the model should the return the specified responses
    % columnwise...
    for rr = 1:Nr
        Rc{rr}.r = Rci(:,rr);
        Rc{rr}.t = Rtype{rr};
        if exist('f','var')
            Rc{rr}.f = f;
        end
    end
end

end



function cost = costSurr(x,S,OPTopts)

% Function to run the Surrogate model S at x to return the cost as
% calculated from the information in OPTopts.
% S can be a cell array of Surrogates, each returning a different response
% at x, with the length of the array equal to the length of the type vector
% in OPTopts.Rtype

if length(S) == 1 && ~iscell(S), S = {S}; end
Nr = length(S); % Number of responses
x = reshape(x,length(x),1);

Rs = cell(1,Nr);
for rr = 1:Nr
    Rs{rr}.r = evalSurr(x,S{rr});
    Rs{rr}.t = OPTopts.Rtype{rr};
    if isfield(S{rr},'f'), Rs{rr}.f = S{rr}.f; end
end
cost = costFunc(Rs,OPTopts);

% if cost == 0, keyboard; end

end

function cost = costFunc(R,GOALS)

% Function to calculate the cost of a given response R as specified by the
% goals in GOALS
% R is a cell array of structures containing the response in R.r, the type R.t, and the
% (optional) domain (typically frequency) in R.f.
% R can also be a structure if only one type of response is considered.
% GOALS can contain (typically a subset of OPTopts used in the main function):
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11dB'
%   goalType:   Cell array of goal types {1,Ng}
%               Valid types:
%                   'lt' (Less than)
%                   'gt' (Greater than) 
%                   'eq' (equal to) - todo
%                   'minimax' 
%   goalVal:    Cell array of goal values {1,Ng} - same order as goalType
%   goalWeight: Vector of goal weights [1,Ng] (default equal weights)
%   goalStart:  Cell array of start of valid goal domain {1,Ng} (optional)
%   goalStop:   Cell array of stop of valid goal domain {1,Ng} (optional)
%   errNorm:    Cell array of type of error norms to use for optimization {1,Ng}
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)
% Todo:
%     - More goalType functions:
%     'gt', 'eq', 'minimax', etc...


% Make R a cell array if only one structure is passed.
if length(R) == 1 && ~iscell(R), R = {R}; end

Nr = length(R);
Ng = length(GOALS.goalType);
[cSum,wSum] = deal(0);

for gg = 1:Ng
    goalType = GOALS.goalType{gg};
    G.goalVal = GOALS.goalVal{gg};
    if isfield(GOALS,'goalStart'), G.goalStart = GOALS.goalStart{gg}; end
    if isfield(GOALS,'goalStop'), G.goalStop = GOALS.goalStop{gg}; end
    G.errNorm = 'L1';
    if isfield(GOALS,'errNorm'), G.errNorm = GOALS.errNorm{gg}; end
    goalWeight = 1;
    if isfield(GOALS,'goalWeight'), goalWeight = GOALS.goalWeight{gg}; end
    wSum = wSum + goalWeight;
    
    % Special case for complex and dB S11 goals...
    if strcmp(R{1}.t,'S11complex') && strcmp(GOALS.goalResType{gg},'S11dB')
        Ri = R{1};
        Ri.r = dB20(R{1}.r);
    else
        Ri = R{1};
    end
    tt = 1;
    while tt < Nr
        if strcmp(R{tt}.t,GOALS.goalResType{gg})
            Ri = R{tt};
            break;
        end
        % Special case for complex and dB S11 goals...
        if strcmp(R{tt}.t,'S11complex') && strcmp(GOALS.goalResType{gg},'S11dB')
            Ri = R{1};
            Ri.r = dB20(R{1}.r);
            break;
        end
        tt = tt + 1;
    end
    if strcmp(goalType,'lt')
        c0 = costlt(Ri,G);
    elseif strcmp(goalType,'gt')
        c0 = costgt(Ri,G);
    elseif strcmp(goalType,'minimax')
        c0 = costminimax(Ri,G);
    end
    cSum = cSum + goalWeight*c0;
end
cost = cSum/wSum;

end

function c0 = costlt(R,G)
% Function to evaluate an 'lt' (less than) goal function for the
% given response R.
% R is a structure containing the response in R.r and the
% (optional) domain (typically frequency) in R.f.  If the domain is
% not specified the (optional) goal limits will be interpreted as indeces in
% the vector.
% The structure G contains:
%   goalVal:    Goal value
%   goalStart:  Start of valid goal domain (optional)
%   goalStop:   Stop of valid goal domain (optional)
%   errNorm:    Type of error norm to use for optimization
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)

Nm = length(R.r);
if ~isfield(R,'f')
    if isfield(G,'goalStart')
        iStart = G.goalStart;
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = G.goalStop;
    else
        iStop = Nm;
    end
else
    if isfield(G,'goalStart')
        iStart = find(R.f >= G.goalStart,1);
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = find(R.f <= G.goalStop,1,'last');
    else
        iStop = Nm;
    end
end

Rvalid = R.r(iStart:iStop);
y = Rvalid - G.goalVal;
y(y < 0) = 0;
c0 = Lnorm(y,G.errNorm);
end

function c0 = costgt(R,G)
% Function to evaluate an 'gt' (greater than) goal function for the
% given response R.
% R is a structure containing the response in R.r and the
% (optional) domain (typically frequency) in R.f.  If the domain is
% not specified the (optional) goal limits will be interpreted as indeces in
% the vector.
% The structure G contains:
%   goalVal:    Goal value
%   goalStart:  Start of valid goal domain (optional)
%   goalStop:   Stop of valid goal domain (optional)
%   errNorm:    Type of error norm to use for optimization
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)

Nm = length(R.r);
if ~isfield(R,'f')
    if isfield(G,'goalStart')
        iStart = G.goalStart;
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = G.goalStop;
    else
        iStop = Nm;
    end
else
    if isfield(G,'goalStart')
        iStart = find(R.f >= G.goalStart,1);
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = find(R.f <= G.goalStart,1,'last');
    else
        iStop = Nm;
    end
end

Rvalid = R.r(iStart:iStop);
y = Rvalid - G.goalVal;
y(y > 0) = 0;
c0 = Lnorm(y,G.errNorm);
end


function c0 = costminimax(R,G)
% Function to evaluate an 'minimax' goal function for the
% given response R.
% R is a structure containing the response in R.r and the
% (optional) domain (typically frequency) in R.f.  If the domain is
% not specified the (optional) goal limits will be interpreted as indeces in
% the vector.
% The structure G contains:
%   goalStart:  Start of valid goal domain (optional)
%   goalStop:   Stop of valid goal domain (optional)
%   errNorm:    Type of error norm to use for optimization
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)

Nm = length(R.r);
if ~isfield(R,'f')
    if isfield(G,'goalStart')
        iStart = G.goalStart;
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = G.goalStop;
    else
        iStop = Nm;
    end
else
    if isfield(G,'goalStart')
        iStart = find(R.f >= G.goalStart,1);
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = find(R.f <= G.goalStop,1,'last');
    else
        iStop = Nm;
    end
end

Rvalid = R.r(iStart:iStop);
c0 = max(Rvalid);
end


function c0 = costeq(R,G)
% Function to evaluate an 'eq' (equal to) goal function for the
% given response R.
% R is a structure containing the response in R.r and the
% (optional) domain (typically frequency) in R.f.  If the domain is
% not specified the (optional) goal limits will be interpreted as indeces in
% the vector.
% The structure G contains:
%   goalVal:    Goal value - can be scalar or vector of the same length as
%               R.r
%   goalStart:  Start of valid goal domain (optional)
%   goalStop:   Stop of valid goal domain (optional)
%   errNorm:    Type of error norm to use for optimization
%               Valid types:
%                   'L1' (L1 norm - default)
%                   'L2' (L2 norm)

Nm = length(R.r);
if ~isfield(R,'f')
    if isfield(G,'goalStart')
        iStart = G.goalStart;
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = G.goalStop;
    else
        iStop = Nm;
    end
else
    if isfield(G,'goalStart')
        iStart = find(R.f >= G.goalStart,1);
    else
        iStart = 1;
    end
    if isfield(G,'goalStop')
        iStop = find(R.f <= G.goalStart,1,'last');
    else
        iStop = Nm;
    end
end

Rvalid = R.r(iStart:iStop);
if length(G.goalVal) == length(R.r)
    y = Rvalid - G.goalVal(iStart:iStop);
else
    y = Rvalid - G.goalVal;
end
c0 = Lnorm(y,G.errNorm);
end


function Ln = Lnorm(y,L)
% Calculates the norm of vector y.  L = 'L1' returns the L1 norm and L
% = 'L2' the L2 norm
if strcmp(L,'L1'), Lp = 1;
elseif strcmp(L,'L2'), Lp = 2; end
Ln = sum(abs(y).^Lp)./length(y);
end


