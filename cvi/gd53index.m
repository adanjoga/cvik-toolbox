function f = gd53index(clust,X,varargin)
% GD53INDEX Evaluation based on the Generalized Dunn-53 criterion.
%   GD53INDEX(CLUST, X) computes the Generalized Dunn-53 criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Generalized Dunn-53 index uses
%   the Euclidean distance between points in X.
%
%   V = GD53INDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Generalized Dunn-53 index.
%   
%   V = GD53INDEX(..., 'DISTANCE', value) computes the Generalized Dunn-53 index using
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
%   eva   = gd53index(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBCVINDEX
%
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
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)'));          
end
den = 2*(sumD./Nk);

% Inter-cluster dispersion
num = inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        DxS = feval(pfun,X(clust==i,:)',M(j,:)');
        DyT = feval(pfun,X(clust==j,:)',M(i,:)');
        num(i,j) = (1/(Nk(i)+Nk(j)))*(sum(DxS)+sum(DyT));
    end
end
% CVI Computation
f = min(num(:)/max(den));