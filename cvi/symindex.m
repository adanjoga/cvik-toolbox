function f = symindex(clust,X,varargin)
% SYMINDEX Evaluation based on the  Symmetry (Sym) criterion.
%   SYMINDEX(CLUST,X) computes the Symmetry criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Symmetry index uses the Symetry
%   distance and internally uses the Euclidean distance.
%
%   V = SYMINDEX(CLUST,X) returns a positive numeric value corresponding to
%   the Symmetry index.
%   
%   V = SYMINDEX(..., 'DISTANCE', value) computes the Symmetry index using
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
%   eva   = symindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   S. Bandyopadhyay and S. Saha, "A Point Symmetry-Based Clustering Technique 
%   for Automatic Evolution of Clusters," 
%   IEEE Transactions on Knowledge and Data Engineering, 
%   Vol. 20, No. 11, pp. 1441–1457, 2008.
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
wcd = zeros(1,K);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    
    dps = symdist(X(clust == i,:),M(i,:),Nk(i),pfun);
    wcd(i) = sum(dps);
end
% Dispersion inter-cluster
Mij = feval(pfun,M',M');
bcd = max(Mij(:));

% CVI Computation
f = bcd/(K*sum(wcd));