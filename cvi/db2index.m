function f = db2index(clust, X,varargin)
% DB2INDEX Evaluation based on the Enhanced Davies-Bouldin criterion.
%   DB2INDEX(CLUST, X) computes the Enhanced Davies-Bouldin criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Enhanced Davies-Bouldin index uses
%   the Euclidean distance between points in X.
%
%   V = DB2INDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Enhanced Davies-Bouldin index.
%   
%   V = DB2INDEX(..., 'DISTANCE', value) computes the Enhanced Davies-Bouldin index using
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
%   eva   = db2index(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   Minho Kim and R. S. Ramakrishna, "New Indices for Cluster Validity Assessment," 
%   Pattern Recognition Letters, 
%   Vol. 26, No. 15, pp. 2353-2363, 2005.
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

% Intra-cluster cohesion (compactness)
dmin = feval(pfun,M',M');
dmin(logical(eye(K))) = Inf;
dmax = -inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        dmax(i,j) = Si(i)+Si(j);
    end
end
% CVI Computation
f = mean(max(dmax,[],1)./min(dmin,[],1));