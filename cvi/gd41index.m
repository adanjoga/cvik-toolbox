function f = gd41index(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -inf;
    return;
end

% Dispersion inter-cluster
M = NaN(K,size(X,2));
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);  
end
num = feval(pfun,M',M');
num(logical(eye(K))) = Inf;

% Dispersion intra-cluster
den = zeros(1,K);
for i = 1:K
    den(i) = max(max(DXX(clrs==i,clrs==i)));
end

% Evaluacion del IVG
f = min(num(:)/max(den));