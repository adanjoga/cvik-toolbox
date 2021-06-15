function [IDX,Mb,bFit,mFit,gfit] = acde(X,Kmax,NP,Gmax,cvi,distance)
% ACDE An automatic clustering algorithm based on differential evolution.
%   IDX = ACDE(X, KMAX, NP, GMAX, CVI) partitions the data points in X into
%   K clusters. ACDE is an evoluationary clustering algorithm which
%   automaticly finds a clustering with K cluster, where K = [2,KMAX].
%   value which can be used for estimating the number of clusters on data.
%
%   X must be an N-by-P matrix of data points with one row per observation 
%   and one column per variable. KMAX is an integer value representing 
%   the maximum number of cluster. NP and GMAX are the population size and  
%   the maximum number of generations, respectively, of the evolutionary 
%   clustering algorithm ACDE. CVI is a string representing the criterion
%   to be used as an objective funtion. See the list of available CVIs 
%   by typing 'help cviconfig'.
%
%   ACDE returns an N-by-1 vector IDX containing the cluster indices
%   of each point. By default, ACDE uses the Euclidean distance.
%
%   [IDX, MB] = ACDE(X, KMAX, NP, GMAX, CVI) returns the K cluster centroid 
%   locations in the K-by-P matrix MB.
%
%   [IDX, MB, BFIT, MFIT] = ACDE(X, KMAX, NP, GMAX, CVI) returns the best
%   and average fitness values, respectively, during the clustering task 
%   in the GMAX-by-1 vectors BFIT and MFIT.
%
%   [IDX, MB, BFIT, MFIT, GFIT] = ACDE(X, KMAX, NP, GMAX, CVI) returns the 
%   global fitness value in the GFIT variable.
%   
%   [...] = ACDE(..., DISTANCE) indicates the distance measure that ACDE 
%   should consider for the computation of the CVI or objective funtion.
%   The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%       'scorr'         - Spearman's correlation coefficient.
%       'lap'           - Laplacian distance.
%
%
%   Example:
%   -------
%   % Parameters related to the automatic clustering problem
%   D = load('Data_4_3');
%   X = D.data(:,1:end-1);
%   Kmax    = 10;               % Maximum number of clusters
%   NP      = 10*size(X,2);     % Population size
%   Gmax    = 100;              % Number of generations
%   CVI     = 'ch';             % CVI name
%   Dist    = 'euc';            % Distance funtion
%
%   % Run the ACDE algorithm
%   [Yb,Pb,bFit] = acde(X, Kmax, NP, Gmax, CVI, Dist);
%
%
%   See also EVALCVI, CVICONFIG
%
%
%   Reference:
%   ----------
%   D. Das, A. Abraham, A. Konar, "﻿Automatic Clustering Using an Improved 
%   Differential Evolution Algorithm," IEEE Transactions on Systems, Man 
%   and Cybernetics, Vol. 38, No. 1, pp. 218-237, 2008.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

% -------------- TO BE REMOVED
% Xmean = mean(X,1)';     % Global mean value of input data X
% DXX = feval(pfun,X',X');
% if strcmpi(cvi,'sdbw')
%     sX  = var(X,0,1);   % Global variance of input data X
%     sXX = sqrt(sum(sX.^2,2)); % 2- norm of input data X
% end
% -------------------------------------------------------------------------

if nargin < 6 || isempty(distance)
    distance = 'euc';
end

[cvifun,opt] = cviconfig(cvi);
pfun = proxconfig(distance);

% Variables of the ACDE algorithm to adjust the crossover probability (CR)
TH = 0.5;
CRmax = 1.0;
CRmin = 0.5;
G = Gmax;
% Generation of the initial population, centroids, and threshold values
D = size(X,2);
Xmin = min(X,[],1);
Xmax = max(X,[],1);
mnmx = Xmin+(Xmax-Xmin);
M = cat(2,repmat(mnmx,[Kmax 1 NP]).*rand(Kmax,D,NP), rand(Kmax,1,NP));

bFit = zeros(1,G);
mFit = zeros(1,G);
fit = zeros(1,NP);

% Evaluation of the initial population
for i = 1:NP
    % Validation of the minimum number of active centroids
    Ti = checkTi(M(:,end,i),Kmax,TH);
    Ki = sum(Ti>TH);
    Mi = M(Ti>TH,1:end-1,i);
    
    % Validation of the minimum number of data points per cluster (>=2)
    [Mi, clrs] = checkMi(X,Mi,Ki,pfun);
    
    % Evaluation of a candidate clustering solution
    fit(i) = feval(cvifun,X,clrs,'distance',distance);

    % Update the corresponding vector in the initial population
    M(Ti>TH,1:end-1,i) = Mi;
    M(:,end,i) = Ti;
end

% Main loop of ACDE algorithm
for g = 1:Gmax
   % Best individual in the population
   if strcmpi(opt,'mx')
       [~,ib] = max(fit);
   elseif strcmpi(opt,'mn')
       [~,ib] = min(fit);
   end 
   % Best fitness in the current population
   bFit(g) = fit(ib);
   % Average fitness in the current population
   mFit(g) = sum(fit(~isinf(fit)&~isnan(fit)))/NP;
   
   disp(['Iteration: ' num2str(g) '| Best: ' num2str(bFit(g)) '| Mean: ' num2str(mFit(g))]);
   
   % Update the crossover probability for the current population
   CR = CRmax-(CRmax-CRmin)*(g/Gmax);
   
   % For each individual in the population
   for i = 1:NP 
       % Target vector
       Mi = M(:,:,i);
       fMi = fit(i);
       
       % MUTATION (DONOR VECTOR: Vi)
       j = setdiff(randperm(NP),i,'stable');
       F = 0.5*(1+rand(1,2));
       Vi = M(:,:,j(1)) + F(1)*(M(:,:,j(2)) - M(:,:,j(3)));
       % CROSSOVER/RECOMBINATION (TRIAL VECTOR: Ri)
       Ri = Mi;
       idx = rand(size(Mi)) < CR;
       Ri(idx) = Vi(idx);
       % Split the individual: centroids and thresholds
       Ti = Ri(:,end);
       Ui = Ri(:,1:end-1);  
       
       % Validate boundary constraints of the centroids
       Ui = boundConstraint(Ui,Kmax,D,Xmin,Xmax);
       % Validate boundary constraints of the thresholds
       Ti = checkTi(Ti,Kmax,TH);
       
       % Get the active centroids
       Ki = sum(Ti>TH); Ut = Ui(Ti>TH,:);
       % Validation of the minimum number of data points per cluster (>=2)
       [Ut,clrs] = checkMi(X,Ut,Ki,pfun);
       % Evaluation of the candidate clustering solution
       fUi = feval(cvifun,X,clrs,'distance',distance);
       % Update the corresponding individual
       Ui(Ti>TH,:) = Ut; Ri = [Ui,Ti];
       
       % SELECTION (elitist replacement)
       if (strcmpi(opt,'mx')&&(fUi >= fMi))||(strcmpi(opt,'mn')&&(fUi <= fMi))
           M(:,:,i) = Ri;
           fit(i) = fUi;
       end
   end
end
% Selection of the best individual
if strcmpi(opt,'mx')
   [~,ib] = max(fit);
elseif strcmpi(opt,'mn')
   [~,ib] = min(fit);
end
bFit(G) = fit(ib);
gfit    = fit(ib);

Tb = M(:,end,ib);
Mb = M(Tb>TH,1:end-1,ib);


if (strcmpi(cvi,'sym')||strcmpi(cvi,'sdunn')||strcmpi(cvi,'sdb'))
    IDX = symclrs(X,Mb,pfun);
else
    D = feval(pfun,X',Mb');
    [~,IDX] = min(D,[],2);
end

% ---------------------------------------------------------------------
function [Mi,clrs] = checkMi(X,Mi,Ki,pfun)
DXM = feval(pfun,X',Mi'); % Coomputation of the distance between points and centroids
[~,clrs] = min(DXM,[],2); % Assignation of each data point to their nearest centroids

N  = numel(clrs);                       % Number of data points
Nk = accumarray(clrs,ones(N,1),[Ki,1]); % Data points assigned to clusters

% Each clusted should contain at least two data points
if sum(Nk<2)
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
    
    % Computation of the new clustering solution based on the new centroids
    DXM = feval(pfun,X',Mi');
    [~,clrs] = min(DXM,[],2);
end

%---------------------------------------------------------------------
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

%---------------------------------------------------------------------
function Y = symclrs(X,Mi,pfun)
K = size(Mi,1);
N = size(X,1);

PSXM = NaN(N,K);
DXM = feval(pfun,X',Mi');
M2 = 2*Mi;
for i = 1:K
    % The symmetrical or reflected point
    rpx = repmat(M2(i,:),N,1)-X; 
    % The knear=2 unique nearest neighbors
    dsym = feval(pfun,rpx',X');
    dsym = sort(dsym,2);
    dsym = mean(dsym(:,1:2),2);
    % computation of PS distance
    PSXM(:,i) = dsym .* DXM(:,i);
end
[~,Y] = min(PSXM,[],2);
