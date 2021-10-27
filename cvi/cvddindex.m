function f = cvddindex(clust,DXX) 
% CVDDINDEX Evaluation based on the CVDD index.
%   CVDDINDEX(CLUST, DXX) computes the CVDD criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = CVDDINDEX(CLUST, DXX) returns a positive numeric value corresponding to
%   the CVDD index.
%   
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = cvddindex(clust, DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, CVNNINDEX, DUNNINDEX, CINDEX
%
%
%   Reference:
%   ----------
%   Lianyu Hu and Caiming Zhong.  
%   "An internal validity index based on density-involved distance,"
%   IEEEAccess, vol 7, pp. 40038–40051, 2019.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
%   Note: The main code was written by  
%   Lianyu Hu, Department of Computer Science, Ningbo University 
%   February 2019
% ------------------------------------------------------------------------

%Validation of the clustering solution
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = -inf;
    return;
end

% WARNING: The DDX could be precomputed before hand
% Computation of the density-involved distance of DXX
DDX = Density_involved_distance(DXX, 7); 

sc_list = zeros(K,1);  % separation
com_list = zeros(K,1); % compactness

for i = 1: K
    a = find(clust == i);
    b = clust ~= i;
    n = length(a);
    if isempty(a)~=1
        %  compute the separation sep[i]
        sc_list(i,1) = min(min(DDX(a,b)));
        %  compute the compactness com[i]
        try
            Ci = fast_PathbasedDist(DXX(a,a));
            com_list(i,1) = (std2(Ci)/n)*mean(Ci(:));
        catch
            com_list(i,1) = max(com_list);
        end
    else
        sc_list(i,1)=0;
        com_list(i,1) = max(com_list);
    end
end
% compute the validity index CVDD 
sep = sum(sc_list);
com = sum(com_list);
try
    f = sep/com;
catch
    % because of 0 distances between some points (total overlap)
    f = 0;
end
end