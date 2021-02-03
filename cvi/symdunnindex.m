function f = symdunnindex(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -Inf;
    return;
end
% Evaluacion del IVG
M = NaN(K,size(X,2));
dmax = zeros(1,K);
dmin = zeros(K);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    
    dps = symdist(X(clrs == i,:),M(i,:),Nk(i),pfun);
    dmax(i) = max(dps);    
    for j = setdiff(1:K,i)
        dmin(i,j) = min(min(DXX(clrs == i,clrs==j)));
    end
end
f = min(dmin(~eye(K))/max(dmax));