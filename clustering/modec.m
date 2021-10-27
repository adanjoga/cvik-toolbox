%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% When one objective is optimized, the standard DE runs; if two or more
% objectives are optimized, the greedy selection step in DE algorithm is 
% performed using a dominance relation.
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341-359.
% 
% Das, S., Suganthan, P. N., 2010. Differential evolution: A survey of the 
% state-of-the-art. IEEE Transactions on Evolutionary Computation. Vol 15,
% 4 - 31.

%% MODE Algoritm
function OUT=modec(MODEData)
%% Reading parameters from MODEData
global X DXX
MAXGEN      = MODEData.MAXGEN;      % Maximum number of generations.
POPSIZE     = MODEData.POPSIZE;     % Population size.
Nobj        = MODEData.NOBJ;        % Number of objectives.

X           = MODEData.X;           % Dataset.
Kmax        = MODEData.KMAX;        % Maximum number of clusters.
CVIs        = MODEData.CVIs;        % Instance.

%% Initialize and evaluate population
pfun    = proxconfig('euc');    % Defines the proximity measure function.
DXX     = real(feval(pfun,X',X'));    % Matrix of distances based on pfun.

[cvifun1,~] = cviconfig(CVIs{1});
[cvifun2,~] = cviconfig(CVIs{2});

% Fixed variables related to the JADE
TH = 0.5; 
c = 0.1;   % Constante positiva entre [0,1]
CRm = 0.5; % Mu CR
Fm  = 0.5; % Mu F

% initialization of the parent population
d = size(X,2);
Xmin = min(X,[],1);
Xmax = max(X,[],1);
mnmx = Xmin+(Xmax-Xmin);
C = repmat(mnmx,[Kmax 1 POPSIZE]).*rand(Kmax,d,POPSIZE);  % Centroides
T = rand(Kmax,1,POPSIZE);                              % Thresholds

Parent = cat(2,C,T);            % Parent population.
FES    = 0;                     % Function Evaluation.

JxParent = zeros(POPSIZE,Nobj); % Fitness of parent population.
CLRs = nan(size(X,1),POPSIZE);  % Clustering solutions
% evaluatation of parent population
for i = 1:POPSIZE
    % Verify the number of active prototypes
    Ti = checkTi(Parent(:,end,i),Kmax,TH);
    Ta = Ti > TH;                      
    Ki = sum(Ta);               % Number of active prototypes
    Pi = Parent(Ta,1:end-1,i);   % Gets the matrix of active prototypes
    
    % Verify the number points (at least 3 points) per group
    [Mi,clrs] = checkMi(X,Pi,Ki,pfun); % Gets Ci = Centroids, Mi = Medoids
    
    f1 = feval(cvifun1,clrs,X);
    f2 = feval(cvifun2,clrs,X);
    JxParent(i,:) = [f1 f2];
    
    Parent(Ta,1:end-1,i) = Mi;	% Update the prototypes
    Parent(:,end,i)      = Ti;	% Update the thresholds
    CLRs(:,i) = clrs;
end
FES = FES+POPSIZE;   

%% Evolution process
for g=1:MAXGEN 
    disp(['Iteration ' num2str(g)]);
    Scr = [];
    Sf = [];
    % Genera CRi y Fi para toda la poblacion
    [Fi, CRi] = randFCR(POPSIZE,CRm,0.1,Fm,0.1);
    for i=1:POPSIZE
        Mi = Parent(:,:,i);
        j = setdiff(randperm(POPSIZE),i,'stable');
                
        % Mutant vector calculation DE/current-to-rand/1
        Vi = Mi + Fi(i)*(Parent(:,:,j(1)) - Mi) + Fi(i)*(Parent(:,:,j(2)) - Parent(:,:,j(3)));
 
        % Crossover operator
        Ui = Mi;
        jrand = randi(d*Kmax,1);
        irand = rand(size(Mi))<CRi(i);
        irand(jrand) = 1;
        Ui(irand) = Vi(irand);
       
        % Split: prototypes and thresholds
        Ui = Ui(:,1:end-1);
        Ti = Ui(:,end);
        
        Ui = boundConstraint(Ui,Kmax,d,Xmin,Xmax); % Bounds of prototypes
        Ti = checkTi(Ti,Kmax,TH);                  % Bounds of thresholds
        
        % Verify the number of active prototypes
        Ta = Ti > TH;
        Ki = sum(Ta);
        Ut = Ui(Ta,:);
        
        % Verify the number points (at least 3 points) per group
        [Ut,clrs] = checkMi(X,Ut,Ki,pfun);
        
        f1 = feval(cvifun1,clrs,X);
        f2 = feval(cvifun2,clrs,X);
        fUi = [f1 f2];
        
        % Update the modified vector and join the thresholds
        Ui(Ta,:) = Ut;
        Ui = cat(2,Ui,Ti);
        
        % Selection in DEMO algorithm
        if fUi <= JxParent(i,:)
            Parent(:,:,i) = Ui;
            JxParent(i,:) = fUi;
            CLRs(:,i) = clrs;
            
            Scr = cat(1,Scr,CRi(i));
            Sf  = cat(1,Sf,Fi(i));
        end
    end
    
    % Update the control paramters mu CR and mu F
    if ~isempty(Scr) && sum(Sf) > 0
        CRm = (1 - c) * CRm + c * mean(Scr);
        Fm  = (1 - c) * Fm + c * sum(Sf .^ 2) / sum(Sf); % Lehmer mean
    end 
    
    FES=FES+POPSIZE;

	PSet    = Parent;
    PFront  = JxParent;

    OUT.XPOP           = Parent;    % Population
    OUT.FPOP           = JxParent;  % Poopulation's Objective Vector
    OUT.PSet           = PSet;      % Pareto Set
    OUT.PFront         = PFront;    % Pareto Front
    OUT.Param          = MODEData;  % MODE Parameters
    OUT.CLRs           = CLRs;      % MODE Clusterings
    MODEData.CounterGEN = g;
    MODEData.CounterFES = FES;
    OUT=PrinterDisplay(OUT,MODEData); % To print results on screen
    
end

% Report solutions
OUT.XPOP = PSet;
OUT.FPOP = PFront;
[OUT.PFront, OUT.PSet, OUT.PClrs] = DominanceFilter(PFront,PSet,CLRs); %A Dominance Filter

% disp('-------------------------------------------------------------------')
% disp('Red  asterisks : Set Calculated.')
% disp('Black circles : Pareto Front.')
% disp('Blus square : Pareto Front.')
% disp('-------------------------------------------------------------------')

F=OUT.PFront;
for i=1:size(F,1)   
    figure(123); hold on;
    plot(F(i,1),F(i,2),'ok','MarkerFaceColor','k'); grid on; hold on;
end
%plot([MODEData.bestValue,MODEData.bestValue],[0,max(OUT.FPOP(:,2))],'--');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% Print and Display information
function [OUT]=PrinterDisplay(OUT,Data)

if mod(Data.CounterGEN,1)==0
    figure(123);
    plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r'); grid on;
end

%% Dominance Filter
function [PFRONT,PSET,PClrs] = DominanceFilter(F,C,CLRs)

NPOP = size(F,1);
PFRONT  = [];
PSET    = [];
PClrs   = [];
k = 0;
for i = 1:NPOP
    Dominado = 0;
    for j = 1:NPOP
        if F(i,:) == F(j,:)
            if i > j
                Dominado = 1;
                break;
            end
        else
            if F(i,:) >= F(j,:)
                Dominado = 1;
                break;
            end
        end
    end
    
    if Dominado == 0
        k = k+1;
        PFRONT = cat(1,PFRONT,F(i,:));
        PSET = cat(3,PSET,C(:,:,i));
        PClrs = cat(2,PClrs, CLRs(:,i));
    end
end

% WARNING!!! This sorting works for 2D fronts in a minimization problem
[~,Idx]=sort(PFRONT(:,1));
PFRONT = PFRONT(Idx,:);
PSET = PSET(:,:,Idx);
PClrs = PClrs(:,Idx);

%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Parameter adaptation function
function [F,CR] = randFCR(NP,CRm,CRsigma,Fm,Fsigma)
CR = CRm + CRsigma * randn(NP, 1);
CR = min(1, max(0, CR));                % truncated to [0 1]
% generate F
F = randCauchy(NP, 1, Fm, Fsigma);
F = min(1, F);                          % truncation
% we don't want F = 0. So, if F<=0, we regenerate F (instead of trucating it to 0)
pos = find(F <= 0);
while ~ isempty(pos)
    F(pos) = randCauchy(length(pos), 1, Fm, Fsigma);
    F = min(1, F);                      % truncation
    pos = find(F <= 0);
end

%---------------------------------------------------------------------
%% Cauchy distribution: cauchypdf = @(x, mu, delta) 1/pi*delta./((x-mu).^2+delta^2)
function result = randCauchy(m, n, mu, delta)
% http://en.wikipedia.org/wiki/Cauchy_distribution
result = mu + delta * tan(pi * (rand(m, n) - 0.5));
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% Validation of active prototypes
function Ti = checkTi(Ti,K,TH)
% Validation lower and upper bounds, values shold be in the range [0,1]
Ti(Ti > 1) = 1;
Ti(Ti < 0) = 0;

% Validation of the minimum number of active thresholds, it whould be two (2)
a = TH; b = 1;
if sum(Ti>TH) < 2
    rp = randperm(K);
    Ti(rp(1)) = a + (b-a).*rand();
    Ti(rp(2)) = a + (b-a).*rand();
end

%% Validation of number of points per group
function [Mi,clrs] = checkMi(X,Mi,Ki,pfun)
DXM = feval(pfun,X',Mi'); % Coomputation of the distance between points and centroids
[~,clrs] = min(DXM,[],2); % Assignation of each data point to their nearest centroids

N  = numel(clrs);                       % Number of data points
Nk = accumarray(clrs,ones(N,1),[Ki,1]); % Data points assigned to clusters

% Each clusted should contain at least two data points
if sum(Nk<1)
    % Assign N/Ki data points to each cluster
    tclrs = zeros(N,1);
    Ni = floor(N/Ki);
    idx = 1; 
    for i=1:Ki
        tclrs(idx:Ni*i) = i;  % assign Ni data points to the i-th cluster
        idx = idx + Ni;
    end
    if mod(N,Ki)
        tclrs(idx:end) = Ki; % assign the remaining data points to Ki
    end
    % Update the centroids
    for i=1:Ki
        members = (tclrs == i);
        Mi(i,:) = mean(X(members,:),1);
    end
    
    clrs = tclrs;
    % Computation of the new clustering solution based on the new centroids
    %DXM = feval(pfun,X',Mi');
    %[~,clrs] = min(DXM,[],2);
end
%%