function [Ri,Si,Pi,Ci,Li] = SMmain(xinit,Sinit,SMopts,Mf,Mc,OPTopts)

% Space Mapping main loop

% Inputs:
% xinit:    Initial values of input parameters [Nn,1]
% Sinit:    Initial SM structure (see buildSurr.m/evalSurr.m)
%           Important to include initial implicit parameters
% SMopts:   SM options (see buildSurr.m/evalSurr.m)
% Mf:       Fine model structure
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'FEKO'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% Mc:       Coarse model structure (can be cell array of structures if more than one type has to be calculated to get all the fine model responses)
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'FEKO'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   Iparams:     Cell array of implicit paramater names - same order as xinit {Nn,1}
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
% OPTopts:  Optimization options structure for external loop
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types (so far):
%               'S11dB'
%               'S11complex'
%               'Gen' - generic case for use with MATLAB models 
%   Ni:         Maximum number of iterations
%   globOpt:    Flag to run PBIL (1 for only first iteration, 2 for all iterations) (default 0)
%   M_PBIL:     Vector of bits for the global search variables (see PBILreal.m)
%   globOptSM:  Flag to run PBIL during the PE process (1 for only first iteration, 2 for all iterations) (default 0)
%   goalResType:Cell array of response names to consider for the different goals {1,Ng}
%               Valid types:
%               'S11dB'
%               'S11complex'
%               'Gen' - general            
%   goalType:   Cell array of goal types {1,Ng}
%               Valid types:
%                   'lt' (Less than)
%                   'gt' (Greater than)
%                   'eq' 
%                   'minimax' 
%                   'bw' (Like 'lt' but also maximizes bandwidth - requires goalCent value)
%   goalVal:    Cell array of goal values {1,Ng} - same order as goalType
%   goalWeight: Vector of goal weights [1,Ng]
%   goalStart:  Cell array of start of valid goal domain {1,Ng}
%   goalStop:   Cell array of stop of valid goal domain {1,Ng}
%   goalCent:   Cell array of centre point of goal domain {1,Ng} (used by the 'bw' goalType) (optional)
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
% Last Modified: 2015-06-26
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
% 2015-03-30: Add FEKO functionality to the fine model function
% 2015-04-04: Add FEKO functionality in the coarse model function
% 2015-05-07: Add limits for x/xp in the cost function
% 2015-06-26: Update fineMod to rather switch between solvers
%             Update coarseMod to rather switch between solvers 
%             Removed costFunc and dependents from function to be used
%             outside as well

% Set defaults
Ni = 10;    % Maximum number of iterations
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
else
    Sinit.xp = [];
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
            % Get the surrogate response after previous iteration
            % optimization - thus at current iteration position
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
    if costFi == 0 && isempty(find(ismember(OPTopts.goalType,'bw'),1))   % Exit if spec is reached (will typically not work for eq and never for minimax, and bw is explicitly excluded)
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
% xi is an array of input parameters - same order as those specified in M
% M is a structure containing all the info to describe to model containing:
%   path:   Full path to file
%   name:   File name of file (without extension)
%   solver:     'CST'/'MATLAB'/'FEKO' (for now)
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

% Call the correct solver
switch M.solver
    
    case 'CST'
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
        
    case 'FEKO'
        % Build parameter string
        parStr = [];
        for nn = 1:Nn
            parStr = [parStr,' -#',M.params{nn},'=',num2str(xi(nn))];
        end
        % Remesh the structure with the new parameters
        FEKOmesh = ['cadfeko_batch ',[M.path,M.name,'.cfx'],parStr];
        system(FEKOmesh)
        % Run FEKO - cannot run with path, so change the directory
        curDir = pwd;
        cd(M.path)
        FEKOrun = ['runfeko ', [M.name,'.cfx']];
        system(FEKOrun)
        cd(curDir)
        % Generate output
        for rr = 1:Nr
            if strncmp(Rtype{rr},'S11',3)
                % Read the S11 touchstone file - must be exported by the FEKO
                % file with the correct name - Name_S11.s1p!
                [Spar,freq] = touchread([M.path,M.name,'_S11.s1p'],1);
                S11 = reshape(Spar(1,1,:),length(freq),1);
                Rf{rr}.f = freq;
            end
            if strcmp(Rtype{rr},'S11dB')
                Rf{rr}.r = dB20(S11);
            elseif strcmp(Rtype{rr},'S11complex')
                Rf{rr}.r = S11;
            end
            Rf{rr}.t = Rtype{rr};
        end
        
    case 'MATLAB'
        Ni = length(M.params);  % This is interpreted as the number of inputs to the function
        inType = [];
        for ii = 1:Ni
            inType = [inType,M.params{ii}];   % Initialise
        end
        switch inType
            case 'xf'
                Rfi = M.name(xi,M.freq);
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
        
    otherwise
        error(['M.solver unknown for fine model evaluation'])
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
%   solver:     'CST'/'FEKO'/'MATLAB' (for now)
%   params:     Cell array of paramater names - same order as xinit {Nn,1}
%   ximin:  Vector of minimum values to limit the parameters [Nn,1]
%   ximax:  Vector of maximum values to limit the parameters [Nn,1]
%   xpmin:  Vector of minimum values to limit the implicit parameters [Nq,1]
%   xpmax:  Vector of maximum values to limit the implicit parameters [Nq,1]
%   freq:       Array of simulation frequencies [Nm,1] (optional)
%   Rtype:      Type of response (cell array if more than one needed)
%               Valid types:
%               'S11dB'
%               'S11complex'
%               'Gen'

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
Nq = length(xp);
% Get number of responses requested
if length(M.Rtype) == 1 && ~iscell(M.Rtype)
    Rtype = {M.Rtype};  % Special case - make a cell
else
    Rtype = M.Rtype;
end
Nr = length(Rtype);
Rc = cell(1,Nr);

% Call the correct solver
switch M.solver
    case 'CST'
        error('CST solver not implimented yet for coarse model evaluations')
        
    case 'FEKO'
        % Build parameter string
        parStr = [];
        for nn = 1:Nn
            parStr = [parStr,' -#',M.params{nn},'=',num2str(xi(nn))];
        end
        % Also include possible implicit parameters
        for qq = 1:Nq
            parStr = [parStr,' -#',M.Iparams{qq},'=',num2str(xp(qq))];
        end
        
        % Remesh the structure with the new parameters
        FEKOmesh = ['cadfeko_batch ',[M.path,M.name,'.cfx'],parStr];
        system(FEKOmesh)
        % Run FEKO - cannot run with path, so change the directory
        curDir = pwd;
        cd(M.path)
        FEKOrun = ['runfeko ', [M.name,'.cfx']];
        system(FEKOrun)
        cd(curDir)
        % Generate output
        for rr = 1:Nr
            if strncmp(Rtype{rr},'S11',3)
                % Read the S11 touchstone file - must be exported by the FEKO
                % file with the correct name - Name_S11.s1p!
                [Spar,freq] = touchread([M.path,M.name,'_S11.s1p'],1);
                S11 = reshape(Spar(1,1,:),length(freq),1);
                Rc{rr}.f = freq;
            end
            if strcmp(Rtype{rr},'S11dB')
                Rc{rr}.r = dB20(S11);
            elseif strcmp(Rtype{rr},'S11complex')
                Rc{rr}.r = S11;
            end
            Rc{rr}.t = Rtype{rr};
        end
        
    case 'AWR'
        error('AWR solver not implimented yet for coarse model evaluations')
        
    case 'ADS'
        error('ADS solver not implimented yet for coarse model evaluations')
        
    case 'MATLAB'
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
    otherwise
        error(['M.solver unknown for coarse model evaluation'])
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

% Check for limits
if isfield(S{1},'ximin')
    ximin = S{1}.ximin;
else
    ximin = -inf.*x;
end
if isfield(S{1},'ximax')
    ximax = S{1}.ximax;
else 
    ximax = inf.*x;
end
% Make dummy xp just for simple limit checking...
if isfield(S{1},'xp') && isfield(S{1},'xpmin') 
    xpmin = S{1}.xpmin;
    xp = S{1}.xp;
else
    xpmin = -inf.*S{1}.xpmin;
    xp = 0;
end
if isfield(S{1},'xp') && isfield(S{1},'xpmax') 
    xpmax = S{1}.xpmax;
    xp = S{1}.xp;
else
    xpmax = -inf.*S{1}.xpmax;
    xp = 0;
end

if all(x <= ximax) && all(x >= ximin) && all(xp <= xpmax) && all(xp >= xpmin)
    Rs = cell(1,Nr);
    for rr = 1:Nr
        Rs{rr}.r = evalSurr(x,S{rr});
        Rs{rr}.t = OPTopts.Rtype{rr};
        if isfield(S{rr},'f'), Rs{rr}.f = S{rr}.f; end
    end
    cost = costFunc(Rs,OPTopts);
else
    cost = inf;
end
% if cost == 0, keyboard; end

end


