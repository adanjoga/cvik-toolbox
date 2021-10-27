function f = lccvindex(clust,X,DXX)
% LCCVINDEX Evaluation based on the Index based on local cores criterion.
%   LCCVINDEX(CLUST,X, DXX) computes the Index based on local cores criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   CLUST is a numeric vector that represents a clustering solution.
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. DXX is an N-by-N dissimilarity matrix.
%
%   V = LCCVINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SDunn index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = LCCVINDEX(CLUST,X, DXX) returns a numeric value corresponding to
%   the index based on local cores.
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = lccvindex(clust, meas, DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, CSINDEX, DUNNINDEX, DBINDEX
%
%   Reference:
%   ----------
%   D. Cheng, Q. Zhu, J. Huang, Q. Wu, and L. Yang.  
%   "A novel clustervalidity index based on local cores."
%   IEEE Transactions on Neural Networks and Learning Systems, 30(4): 985–999, 2019.
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
%   Note: The main code was written by  
%   Dongdong Cheng, Department of Computer Science, Chongqing University 
%   December 2016
% ------------------------------------------------------------------------

% Parameter validation
if nargin < 3
    DXX = pdist2(X,X,'Euclidean'); % Euclidean distance matrix
end
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

% Compute the local cores information (these variables could be pre-computed)
[cores,short_path,local_core] = compute_local_cores(X,DXX);

[~,ncores] = size(cores);
[N,dim] = size(X);
D = zeros(ncores,dim);
cl_cores=zeros(ncores,1);
for i=1:ncores
    D(i,:)=X(cores(i),:);
    cl_cores(i)=clust(cores(i));
end

% Count the number of points belonging to each cores
nl = zeros(ncores,1);
for i=1:ncores
    for j=1:N
        if local_core(j)==cores(i)&&clust(j)>0
            nl(i)=nl(i)+1;
        end
    end
end
% Compute the Silhouette value for each core cluster
[~,s] = computeSWC(D,cl_cores,short_path);
mcv=0;
for i=1:ncores
    mcv=mcv+s(i)*nl(i);
end

f=mcv/N;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% X: The dataset
% clust: The clustering solution
% dist: The graph-based distance between local cores
% -------------------------------------------------------------------------
% Output:
% sil: result the silhoutte width criterion (when using silhouette index)
% s_lccv: the silhouette value of each point (when computing LCCV index)
% -------------------------------------------------------------------------
% Reference: Rousseeuw P. Silhouettes: A graphical aid to the interpretation
%           and validation of cluster analysis. Journal of Computational
%           & Applied Mathematics, 1987, 20(20):53-65.
% -------------------------------------------------------------------------
% Written by Dongdong Cheng
% Department of Computer Science, Chongqing University 
% December 2016
function [sil,s_lccv] =computeSWC(X,clust,dist)

K = numel(unique(clust));

[N,~]=size(X);
cdata=cell(1,K); % the number of points in each cluster
numc=K;
for i=1:K
    nump=0;
    for j=1:N
        if clust(j)==i
            nump=nump+1;
            cdata{1,i}(nump,:)=X(j,:);
            cindex(i,nump)=j;
        end
    end
end

% Don't compute the swc of outliers
numo=0;
if min(clust)<=0
    for i=1:N
        if clust(i)<=0
            numo=numo+1;
        end
    end
end

sil=0;
s_lccv=zeros(N,1);
for i=1:numc
    a=[];
    b=[];
    s=[];
    [np,~]=size(cdata{1,i});
    if np>1
    for j=1:np
        suma=0;
        for k=1:np
            if j~=k
               suma=suma+dist(cindex(i,j),cindex(i,k));
            end
        end
        a(j)=suma/(np-1);
        d=ones(1,numc)*inf;
        for k=1:numc
            if k~=i
                [np2,~]=size(cdata{1,k});
                sumd=0;
                for l=1:np2
                    sumd=sumd+dist(cindex(i,j),cindex(k,l));
                end
                d(k)=sumd/np2;
            end
        end
        b(j)=min(d);
        s(j)=(b(j)-a(j))/max(a(j),b(j));
        s_lccv(cindex(i,j))=s(j);
        sil=sil+s(j); 
    end 
    end
end
sil=sil/(N-numo);
end
