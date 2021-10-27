function f = chindex(clust, X, varargin)
% CHINDEX Evaluation based on the Calinski-Harabasz criterion.
%   CHINDEX(CLUST, X) computes the Calinski-Harabasz criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Calinski-Harabasz index uses
%   the Euclidean distance between points in X.
%
%   V = CHINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Calinski-Harabasz index.
%   
%   V = CHINDEX(..., 'DISTANCE', value) computes the Calinski-Harabasz index using
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
%   eva   = chindex(clust,meas);
%
%
%   See also EVALCVI, CVICONFIG, DBINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%   Reference:
%   ----------
%   T. Calinski, J. Harabasz, "A dendrite method for cluster analysis," 
%   Communications in Statistics - Theory and Methods, 
%   Vol. 3, No. 1, pp. 1-27, 1974.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[mydist] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
pfun = proxconfig(mydist);

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
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);
end
SSW = sum(sumD,1);

% Inter-cluster dispersion
Xmean = mean(X,1)';   % Global mean in X
SSB = sum(Nk .* ((feval(pfun,M',Xmean)).^2));

% CVI Computation
f = ((N-K)/(K-1))*(SSB/SSW);
