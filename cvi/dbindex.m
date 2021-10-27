function f = dbindex(clust,X,varargin)
% DBINDEX Evaluation based on the Davies-Bouldin criterion.
%   DBINDEX(CLUST, X) computes the Davies-Bouldin criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Davies-Bouldin index uses
%   the Euclidean distance between points in X.
%
%   V = DBINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Davies-Bouldin index.
%   
%   V = DBINDEX(..., 'DISTANCE', value) computes the Davies-Bouldin index using
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
%   eva   = dbindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   David L. Davies and Donald W. Bouldin, "A Cluster Separation Measure," 
%   IEEE Transactions on PatternAnalysis and Machine Intelligence, 
%   Vol. 2, No. 1, pp. 224-227, 1979.
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

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)'));          
end
Si = sumD./Nk;

% Inter-cluster dispersion
Mij = feval(pfun,M',M');
R = -inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        R(i,j) = (Si(i)+Si(j))/Mij(i,j);
    end
end

% CVI Computation
f = mean(max(R,[],2));
