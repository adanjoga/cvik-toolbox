function f = xbindex(clust,X,varargin)
% XBINDEX Evaluation based on the Xie-Beni criterion.
%   XBINDEX(CLUST, X) computes the Xie-Beni criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P matrix of data with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Xie-Beni index uses
%   the Euclidean distance between points in X.
%
%   V = XBINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Xie-Beni index.
%   
%   V = XBINDEX(..., 'DISTANCE', value) computes the Xie-Beni index using
%   a specified distance measure. The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%       'scorr'         - Spearman's correlation coefficient.
%       'lap'           - Laplacian distance.
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   eva   = xbindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, DBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   Xuanli Lisa Xie and Gerardo Beni, "A Validity Measure for Fuzzy Clustering," 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%   Vol. 13, No. 8, pp. 841–847, 1991.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

%Parameter validation
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
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);          
end
num = sum(sumD);

% Inter-cluster dispersion
den = feval(pfun,M',M').^2;
den(logical(eye(K))) = Inf;
den = min(den(:));

% CVI Computation
f = (1/N)*(num/den);

