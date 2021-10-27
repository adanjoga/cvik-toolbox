function f = sdbwindex(clust, X,varargin)
% SDBWINDEX Evaluation based on the SDBw criterion.
%   SDBWINDEX(CLUST, X) computes the SDBw criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the SDBw index uses
%   the Euclidean distance between points in X.
%
%   V = SDBWINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SDBw index.
%   
%   V = SDBWINDEX(..., 'DISTANCE', value) computes the SDBw index using
%   a specified distance measure. The available built-in measures are:
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
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   eva   = sdbwindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   M. Halkidi and M. Vazirgiannis, "Clustering Validity Assessment: 
%   Finding the Optimal Partitioning ofa Data Set," 
%   In Proceedings of the IEEE International Conference on Data Mining, 
%   California, USA, pp. 187–194, 2001.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

% Parameter validation
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[Dtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
pfun = proxconfig(Dtype);
% ------------------------------------------------------------------------
% Validation of the clustering solution
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = inf;
    return;
end

% WARNING: If possible compute the sXX variable before calling the SDBw index
% global sXX;
sX  = var(X,0,1);   % Global variance of input data X
sXX = sqrt(sum(sX.^2,2)); % 2- norm of input data X

% Intra-cluster variance
sM  = zeros(K,size(X,2)); % Varianza del cluster Ci
for i = 1:K
    sM(i,:) = var(X(clust==i,:),0,1);
end
sMM = sqrt(sum(sM.^2,2)); % 2-norma |x| = (x.x)^0.5
% Varianza del dataset: sX  = var(X,0,1);
% 2-norma del dataset:  sXX = sqrt(sum(sX.^2,2));
Scat = (1/K)*sum(sMM/sXX);

% Average standard deviation of clusters
stdev = (1/K)*sqrt(sum(sMM));

% Inter-cluster density
M = NaN(K,size(X,2));
Dmin = zeros(N,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    Dmin(members) = feval(pfun,X(members,:)',M(i,:)');
end

aux1 = zeros(1,K);
for i = 1:K
   c = 0;
   Di = density(Dmin(clust==i,:),[],stdev,pfun);
   aux2 = zeros(1,K-1);
   for j = setdiff(1:K,i)
       c = c+1;
       Uij = 0.5*(M(i,:)+M(j,:));
       Dij = density(X(clust==i|clust==j,:),Uij,stdev,pfun);
       Dj  = density(Dmin(clust==j,:),[],stdev,pfun);
       aux2(c) = Dij/max(max(Di,Dj),1);
   end
   %aux2(isinf(aux2)|isnan(aux2)) = 0;
   aux1(i) = sum(aux2);
end
Dens_bw = (1/(K*(K-1)))*sum(aux1);

% CVI Computation
f = Scat + Dens_bw;

if isnan(f)
    inf;
end
%***********************************************************************
function D = density(Z,U,sd,pfun)
if isempty(U)
    D = sum(Z <= sd,1);
else
    D = sum(feval(pfun,Z',U') <= sd,1);
end