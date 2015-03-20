function [c,fval,exitFlag,output] = PBIL(funfcn,M,options)

% function [c,fval,exitFlag,output] = PBIL(funfcn,M,options)
% performs a PBIL optimization on the function @funfcn which must accept a
% chromosome of M bits as the first argument.
%
% 'options' is a structure containing the following optional fields with 
% [defaults] also shown:
% display = 'off';      % Level of display [ {off} | iter | final ]
% algorithm = 'search';   % Specific variation of PBIL to use [{search} | mutate]
% Pop = 50;      % Population size {int > 1}
% alp = 0.1;      % Learning rate {0,1}
% Nbest = 30;     % Number of iterations to continue searching for new best {int > 1}
% Niter = 500;    % Maximum number of iterations {int > 1}
% Nfeval = 10000;  % Maximum number of function evaluations {int > 1}
% % For the mutate algortihm
% mu = 1;         % Number of best solutions {int > 0}
% bet = 0.0075;    % Negative learning rate {0,1}
% nu = 0;         % Number of worst solutions {int >= 0}
% Mp = 0.02;      % Mutation probability {0,1}
% Mbet = 0.05;    % Mutation shift {0,1}
%
% The best chromosome is returned in c, and the function value in fval.
% Optionally the exitFlag is also returned with:
% exitFlag = 1 - ['No improvement for ',num2str(bestCount),' iterations']
% exitFlag = 2 - 'Maximum number of iterations reached'
% exitFlag = 3 - 'Maximum number of function evaluations reached'
% exitFlag = 4 - 'Probability vector converged' (Only for 'search' algorithm)
%
% 'output' is a structure containing the iteration sampled information of
% the search:
% output.P = Pout;                      % Probability vector at each iteration
% output.popSizeUnique = popSizeUnique; % Population size for each iteration
% output.evalNum = evalNum;             % Number of function evaluations performed during each iteration
% output.CFbest = CFbest;               % Best chromosome and fuction value at each iteration
% output.bestCount = bestCount;         % Number of iterations with unchanged best value
% output.iterCount = iterCount;         % Number of iterations
% output.fevalCount = fevalCount;       % Number of function evaluations
% output.outMessage = outMessage;       % Output message


% Date created: 2014-09-??
% Dirk de Villiers
% Last Modified: 2014-09-??
% Updates:



% Set defaults
% Display control
if isfield(options,'display')
    if strcmp(options.display,'off')
        d = 0;
    elseif strcmp(options.display,'iter')
        d = 2;
    elseif strcmp(options.display,'final')
        d = 1;
    else
        d = 0;
    end
else
    d = 0;
end
% Specific variation of PBIL to use {[mutate],search}
if isfield(options,'algorithm'), algorithm = options.algorithm; else algorithm = 'search'; end
% Population size
if isfield(options,'Pop'), Pop = options.Pop; else Pop = 50; end
% Learning rate
if isfield(options,'alp'), alp = options.alp; else alp = 0.1; end
% Number of iterations to continue searching for new best
if isfield(options,'Nbest'), Nbest = options.Nbest; else Nbest = 30; end  
% Maximum number of iterations
if isfield(options,'Niter'), Niter = options.Niter; else Niter = 500; end
% Maximum number of function evaluations
if isfield(options,'Nfeval'), Nfeval = options.Nfeval; else Nfeval = 10000; end

% Set control parameters and check algorithm type input
if strcmp(algorithm,'mutate')
    a = 1;
elseif strcmp(algorithm,'search')
    a = 2;
    s = 1 - (2/sqrt(Pop))^(1/M);    % Search rate
    f = 2*s*alp/(1-2*s*(1-alp));    % f
    t = s + (0.5 - s)/10;           % t for convergence of probability vector
else
    a = 1;
    warning('Unknown algorithm requested - using default: mutate');
    algorithm = 'mutate';
end

% Mutation algorithm specific parameters
if a == 1
    % Number of best solutions
    if isfield(options,'mu'), mu = options.mu; else mu = 1; end
    % Negative learning rate
    if isfield(options,'bet'), bet = options.bet; else bet = 0.0075; end
    % Number of worst solutions
    if isfield(options,'nu'), nu = options.nu; else nu = 1; end
    % Mutation probability
    if isfield(options,'Mp'), Mp = options.Mp; else Mp = 0.02; end
    % Mutation shift
    if isfield(options,'Mbet'), Mbet = options.Mbet; else Mbet = 0.05; end
end


% Initial probability and 'best' chromosome
P = 0.5*ones(1,M);
Ci = zeros(1,M);
for mm = 1:M
    r = rand;
    if r >= P(mm), Ci(mm) = 0; else Ci(mm) = 1; end
end

% Build population and evaluate
C(1,:) = Ci;

CF = [];   % Total population of all chromosomes ever evaluated (columns 1:M) and the solutions in the rest of the columns (M+1)
Cbest = Ci;

header = 'Iteration   Func-count     min f(x)      Best-count';
exitFlag = 0;
iterCount = 0;
bestCount = 0;
if d > 1
    disp(header)
end
while exitFlag == 0
    iterCount = iterCount+1;
    Pout(iterCount,:) = P;
    for pp = 2:Pop
        for mm = 1:M
            r = rand;
            if r >= P(mm), C(pp,mm) = 0; else C(pp,mm) = 1; end
        end
    end
    % Get the unique chromosomes in the population, and the locations of them
    [Cu,locU] = unique(C,'rows');
    [popSizeUnique(iterCount),dummy] = size(Cu); % Store the unique population size at each iteration
    % Positions of Cu to include in function evaluations
    posCalc = ones(popSizeUnique(iterCount),1);
    % Find the values that have been previously computed
    if iterCount > 1,
        [tCF,locCu] = ismember(CF(:,1:M), Cu, 'rows');
        posCalc(locCu(tCF)) = 0;  % Only set positions where previous results are known to 0
    end
    
    % Evaluate the function for all new chromosomes
    evalNum(iterCount) = sum(posCalc);
    Fu = zeros(popSizeUnique(iterCount),1);
    for uu = 1:popSizeUnique(iterCount)
        if posCalc(uu)
            Fu(uu,1) = funfcn(C(locU(uu),:));
%             % Evaluate the function here - example gets the decimal value
%             % of the chromosome
%             sCuW = num2str(C(locU(uu),:));
%             sCu = regexprep(sCuW,'\W','');
%             Fu(uu,1) = 1*bin2dec(sCu);
        else
            % Include previous results
            locCF = min(find(locCu == uu)); % Index of the full set equal to the current index of the iteration unique vector
            Fu(uu,1) = CF(locCF,M+1);
        end
    end
    
    % Build the full set
    CF = [CF;Cu,Fu];
    % Reduce to the unique values
    CF = unique(CF,'rows');
    clear Fu
    % Sort solutions from best to worst (low to high)
    CF = sortrows(CF,M+1);
    CFbest(iterCount,:) = CF(1,:);
    if iterCount > 1
        if CFbest(iterCount,M+1) == CFbest(iterCount-1,M+1)
            bestCount = bestCount + 1;
        else
            bestCount = 0;
        end
    end
    
    % Update P
    for mm = 1:M
        if a == 1
            % Update probability towards best solutions
            P(mm) = P(mm).*(1-alp) + alp./mu.*sum(CF(1:mu,mm));
            
            % Update probability away from worst solutions
            if nu > 0
                [iWn,dummy] = size(CF);
                iW1 = iWn-nu+1;
                P(mm) = P(mm).*(1-bet) + bet./nu.*sum(~CF(iW1:iWn,mm));
            end
            
            % Mutate
            rm = rand;
            rs = rand;
            if rm < Mp
                if rs < 0.5
                    P(mm) = P(mm)*(1-Mbet);
                else
                    P(mm) = P(mm)*(1-Mbet) + Mbet;
                end
            end
        elseif a == 2
            P(mm) = ((1-alp).*P(mm) + alp.*CF(1,mm)).*(1-f) + f/2;
        end
    end
    
    fevalCount = sum(evalNum);
    if bestCount >= Nbest
        exitFlag = 1;
        outMessage = ['No improvement for ',num2str(bestCount),' iterations'];
    end
    if iterCount >= Niter
        exitFlag = 2;
        outMessage = 'Maximum number of iterations reached';
    end
    if fevalCount >= Nfeval
        exitFlag = 3;
        outMessage = 'Maximum number of function evaluations reached';
    end
    if a == 2
        if max(0.5 - abs(P - 0.5)) < t
            exitFlag = 4;
            outMessage = 'Probability vector converged';
        end
    end
    
    if d > 1
        disp(sprintf(' %5.0f        %5.0f     %12.6g         %5.0f', iterCount, fevalCount, CF(1,M+1), bestCount));
    end
end

if d > 0
    disp(' ')
    disp(outMessage)
    disp(' ')
    disp(header)
    disp(sprintf(' %5.0f        %5.0f     %12.6g         %5.0f', iterCount, fevalCount, CF(1,M+1), bestCount));
end

c = CF(1,1:M);
fval = CF(1,M+1);
output.P = Pout;                      % Probability vector at each iteration
output.popSizeUnique = popSizeUnique; % Population size for each iteration
output.evalNum = evalNum;             % Number of function evaluations performed during each iteration
output.CFbest = CFbest;               % Best chromosome and fuction value at each iteration
output.bestCount = bestCount;         % Number of iterations with unchanged best value
output.iterCount = iterCount;         % Number of iterations
output.fevalCount = fevalCount;       % Number of function evaluations
output.outMessage = outMessage;       % Output message



