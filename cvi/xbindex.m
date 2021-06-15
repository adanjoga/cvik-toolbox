function f = xbindex(X,clrs,varargin)
% XBINDEX Evaluation based on the Xie-Beni criterion.
%   XBINDEX(X, CLUST) computes the Xie-Beni criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   X is an N-by-P matrix of data with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Xie-Beni index uses
%   the Euclidean distance between points in X.
%
%   V = XBINDEX(X, CLUST) returns a positive numeric value corresponding to
%   the Calinski-Harabasz index.
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
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   eva   = xbindex(meas, clust);
%
%
%   See also EVALCVI, CVICONFIG, SILINDEX, DUNNINDEX
%
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

% Parameters validations
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

pnames = {'distance', 'datatype'}; pdvals = {'euc', 'feature'};
[Dtype, Xtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});

pfun = proxconfig(Dtype);
if strcmp(Xtype,'relational')
    error('The DataType value for this index must be "feature".')
end

% Validation of the clustering solution 
N = numel(clrs);
clusts = unique(clrs);
K = numel(clusts);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = inf;
    return;
end

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
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

