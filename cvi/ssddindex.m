function f = ssddindex(clust,X,DXX)
% SSDDINDEX Evaluation based on the SSDD criterion.
%   SSDDINDEX(CLUST, X, DXX) computes the SSDD criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   CLUST is a numeric vector that represents a clustering solution.
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. DXX is an N-by-N dissimilarity matrix.
%
%   V = SSDDINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SDunn index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = SSDDINDEX(CLUST, X, DXX) returns a positive numeric value corresponding to
%   the SSDD index.
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = ssddindex(clust, meas, DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, CSINDEX, DUNNINDEX, DBINDEX
%
%   Reference:
%   ----------
%   Shaoyi Liang, Deqiang Han, and Yi Yang. 
%   "Cluster validity index for irregular clustering results."
%   Applied Soft Computing, 95: 106583, 2020
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

% Parameter validation
if nargin < 3
    DXX = pdist2(X,X,'Euclidean'); % Euclidean distance matrix
end
%WARNING1: Each cluster should contain at least 3 backbone points
%WARNING2: 'l' should be approximately 3% of the considered cluster
%WARNING3: the minimum number of points per cluster should be 3 points
% ------------------------------------------------------------------------
%Validation of the clustering solution 
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<3) || K==1
    f = inf;
    return;
end

% CVI computation
d = size(X,2);
DC = zeros(1,K);
ci = cell(1,K);
ICD = zeros(1,K);
BAVDensCi = zeros(1,K);
for k = 1:K
    id = clust==k;
    l = round(0.02*sum(id))+2; % (because of WARNING2, it should be l>1)
    DXK = DXX(id,id);
    Xk = X(id,:);
    % find ci's backbone points 
    % (to consider WARNING1 here!!!)
    [dist,idx] = sort(DXK,2);
    dist(:,1) = []; idx(:,1) = [];
    LD = 1./(mean(dist(:,1:l),2)+eps);
    idx = idx(:,1:l);
    bb = find(sum(LD>LD(idx),2)==l);
    ci{k} = Xk(bb,:);
    % build the MST based on the found backbone points
    D = DXK(bb,bb);
    D = sparse(D); 
    G = graph(D);
    T = minspantree(G);
    list = table2array(T.Edges(:,1));
    diam = table2array(T.Edges(:,2));
    % calculate the BAVDens for all pairs of adjacent vertices in the MST
    BAVDens = zeros(1,size(list,1));
    for i = 1:size(list,1)
        id1 = bb(list(i,1));
        id2 = bb(list(i,2));
        ctr = 0.5*(Xk(id1,:)+Xk(id2,:));
        inside = sum((0.5*diam(i)-pdist2(Xk,ctr,'euclidean'))>=0);
        BAVDens(i) = inside/(diam(i)^d);
    end
    % calculate the density changes DC(ci)
    mx = max(BAVDens);
    mn = min(BAVDens);
    DC(k) = (mx-mn)/mx; 
    BAVDensCi(k) = mean(BAVDens);
    % Density of Region Between One Cluster and All Other Clusters
    idx = cat(2,bb,idx(bb,:));
    ICDk = zeros(1,size(idx,1));
    for i = 1:size(idx,1)
        knnbl = idx(i,:);
        D = DXK(knnbl,knnbl);
        D = sparse(D); 
        G = graph(D);
        T = minspantree(G);
        ICDk(i) = mean(table2array(T.Edges(:,2)));
    end
    ICD(k) = mean(ICDk);
end

flag = zeros(K^2-K,2);
BNPDenspapb = cell(K^2-K,1);
c = 0;
for i = 1:K
   bi = ci{i};
   for j = 1:K
       if i~=j
           c = c+1;
           bj = ci{j};
           D = pdist2(bi,bj,'euclidean');
           [diam,idx] = sort(D,2);
           idx = idx(:,1); diam = diam(:,1);
           bj = bj(idx,:);
           ctr = 0.5*(bi + bj);
           Rij = 0.5*min(ICD(i),ICD(j));
           Dij = (2*Rij)^d;
           BNPDens = zeros(size(ctr,1),1);
           % find all the inter-cluster nearest data pairs between ci and
           % the other clusters
           for k = 1:size(ctr,1)
               inside = (0.5*diam(k)-pdist2(X,ctr(k,:),'euclidean'))>=0;
               Lk = clust(inside);
               cls = unique(Lk)';
               if numel(cls)==2
                  if all(cls==sort([i j]))
                      Xk = X(inside,:);
                      Xi = Xk(Lk==i,:);
                      Xj = Xk(Lk==j,:);
                      di = pdist2(bj(k,:),Xi,'euclidean');
                      [~,idxi] = min(di);
                      pa = Xi(idxi,:);
                      dj = pdist2(bi(k,:),Xj,'euclidean');
                      [~,idxj] = min(dj);
                      pb = Xj(idxj,:);
                      % calculate the BNPDens between all pairs of inter-cluster 
                      % nearest data pairs        
                      ctrk = 0.5*(pa+pb);
                      insidek = (Rij-pdist2(X,ctrk,'euclidean'))>=0;
                      BNPDens(k) = sum(insidek)/Dij;
                  end
               end
           end
           BNPDenspapb{c} = BNPDens;
           flag(c,:) = [i j];

       end
   end
end
% calculate ci’s inner-cluster density and inter-cluster density ratio DR(ci)
DR = zeros(1,K);
for i = 1:K
    % calculate the BNPDens between all pairs of inter-cluster nearest
    % data pairs
    BNPDens = mean(cat(1,BNPDenspapb{flag(:,1)==i}));
    DR(i) = BNPDens/max(BNPDens,BAVDensCi(i));
end
% if nan values, it means that it is not an amogeneous cluster and should
% be penalized by assigning their values to 1 (ssdd is a minimization index)
DR(isnan(DR)) = 1; DC(isnan(DC)) = 1;

% calculate SSDD(C) based on DC(ci) and DR(ci)
alpha = 0.5;
f = mean(alpha*DC + (1-alpha)*DR);
