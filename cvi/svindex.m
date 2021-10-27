function f = svindex(clust,X,varargin)
% SVINDEX Evaluation based on the SV criterion.
%   SVINDEX(CLUST, X) computes the SV criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the SV index uses
%   the Euclidean distance between points in X.
%
%   V = SVINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SV index.
%   
%   V = SVINDEX(..., 'DISTANCE', value) computes the SV index using
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
%   eva   = svindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   K. Zalik and B. Zalik, "Validity Index for Clusters of Different Sizes and Densities," 
%   Pattern Recognition Letters, 
%   Vol. 32, No. 2, pp. 221–234, 2011.
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
    f = -inf;
    return;
end

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
maxD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    maxD(i) = max(feval(pfun,X(members,:)',M(i,:)'));
end
den = sum(maxD,1);

% Inter-cluster dispersion
Cij = feval(pfun,M',M');
Cij(logical(eye(K))) = Inf;
num = sum(min(Cij,[],1));

% CVI Computation
f = num/den;

