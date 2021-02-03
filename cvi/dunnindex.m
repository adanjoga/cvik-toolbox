function f = dunnindex(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -inf;
    return;
end
% Evaluacion del IVG
dmin = inf(K);
dmax = zeros(1,K);
for i = 1:K
    dmax(i) = max(max(DXX(clrs==i,clrs==i)));
    for j = setdiff(1:K,i)
        dmin(i,j) = min(min(DXX(clrs==i,clrs==j)));
    end
end
f = min(dmin(:)/max(dmax));