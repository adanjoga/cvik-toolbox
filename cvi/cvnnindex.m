function f = cvnnindex(clust, DXX)
% CVNNINDEX Evaluation based on the CVNN index.
%   CVNNINDEX(CLUST, DXX) computes the CVNN criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = CVNNINDEX(CLUST, DXX) returns a positive numeric value corresponding to
%   the CVNN index.
%   
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = cvnnindex(clust, DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, CVDDINDEX, DUNNINDEX, CINDEX
%
%
%   Reference:
%   ----------
%   Y. Liu, Z. Li, H. Xiong, X. Gao, J. Wu, and S. Wu, 
%   "Understanding and enhancement of internal clustering validation measures," 
%   IEEE Trans. Cybern., vol. 43, no. 3, pp. 982–994, Jun. 2013.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

%Validation of the clustering solution
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = inf;
    return;
end

% CVI Computation
Knn = 8;
[KNNG]=KNearestNeighborGraph(DXX,Knn);
cluster_weight = zeros(K,1);
compactness = 0;
for i=1:K
    a = find(clust == i);
    b = find(clust ~= i);
    num_a = length(a);
    for n=1:num_a
        if ~isempty(intersect(KNNG{a(n),1},b))
            q = length(intersect(KNNG{a(n),1},b));
            cluster_weight(i) = cluster_weight(i) + q/Knn;
        end     
    end
    cluster_weight(i) = cluster_weight(i)/num_a;
    compactness = compactness + (2/(num_a*(num_a-1)))*sum(sum(DXX(a,a)));
end
separation = max(cluster_weight);
f = separation + compactness;

end