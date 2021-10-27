function f = wbindex(clust,X,varargin)
% WBINDEX Evaluation based on the WB criterion.
%   WBINDEX(CLUST, X) computes the WB criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Davies-Bouldin index uses
%   the Euclidean distance between points in X.
%
%   V = WBINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the WB index.
%   
%   V = WBINDEX(..., 'DISTANCE', value) computes the WB index using
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
%   eva   = wbindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   Q. Zhao and P. Fränti, 
%   "WB-index: A sum-of-squares based index for cluster validity," 
%   Data Knowl. Eng., vol. 92, pp. 77–89, Jul. 2014.
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

% CVI Computation
C_x = mean(X);
W = 0;
B = 0;
for i = 1: K
    C_i = find(clust==i);
    n = length(C_i);
    X_i = X(C_i,:);
    Diff_w = X_i - repmat(mean(X_i,1), length(C_i), 1);
    W = W + sum(sum(Diff_w .* Diff_w, 2));     
    Diff_b = mean(X_i,1)-C_x;
    B = B + n*(sum(Diff_b .* Diff_b));    
end
f = (K*W)/B;
end