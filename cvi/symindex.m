function f = symindex(X,clrs,K,pfun)
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -inf;
    return;
end

% Dispersion intra-cluster (within-cluster dist)
M = NaN(K,size(X,2));
wcd = zeros(1,K);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    
    dps = symdist(X(clrs == i,:),M(i,:),Nk(i),pfun);
    wcd(i) = sum(dps);
end
% Dispersion inter-cluster (between-cluster dist)
Mij = feval(pfun,M',M');
bcd = max(Mij(:));

% Evaluacion del IVG
f = bcd/(K*sum(wcd));