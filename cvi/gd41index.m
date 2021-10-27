function f = gd41index(clust,X,DXX,varargin)
% GD41INDEX Evaluation based on the Generalized Dunn-41 criterion.
%   GD41INDEX(CLUST, X, DXX) computes the Generalized Dunn-41 criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   CLUST is a numeric vector that represents a clustering solution.
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. DXX is an N-by-N dissimilarity matrix.
%
%   V = GD41INDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Generalized Dunn-41 index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = GD41INDEX(CLUST, X, DXX) returns a positive numeric value corresponding to
%   the Generalized Dunn-41 index.
%
%   V = GD33INDEX(..., 'DISTANCE', value) computes the Generalized Dunn-41 index using
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
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = gd41index(clust,meas,DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, CSINDEX, DUNNINDEX, DBINDEX
%
%   Reference:
%   ----------
%   James C. Bezdek and N. R. Pal, "Some New Indexes of Cluster Validity," 
%   IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics), 
%   Vol. 28, No. 3, pp. 301–315, 1998.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
%Parameter validation
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[Dtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
pfun = proxconfig(Dtype);

if nargin < 3
    DXX = real(feval(pfun,X',X'));
end
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

% Dispersion inter-cluster
M = NaN(K,size(X,2));
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);  
end
num = feval(pfun,M',M');
num(logical(eye(K))) = Inf;

% Dispersion intra-cluster
den = zeros(1,K);
for i = 1:K
    den(i) = max(max(DXX(clust==i,clust==i)));
end

% CVI Computation
f = min(num(:)/max(den));