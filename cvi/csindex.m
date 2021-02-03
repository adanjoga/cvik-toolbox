function f = csindex(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);
if numel(clusts) ~= K || sum(Nk<2)
    f = inf;
    return;
end

% Dispersion intra-cluster
M = NaN(K,size(X,2));
num = 0;
for i = 1:K
   XC = DXX(clrs==i,clrs==i);
   num = num+mean(max(XC,[],2));
   
   members = (clrs == clusts(i));
   M(i,:) = mean(X(members,:),1);
end

% Dispersion inter-cluster
Mij = feval(pfun,M',M');
Mij(logical(eye(K))) = Inf;

% Evaluacion del IVG
f = num/sum(min(Mij,[],2));