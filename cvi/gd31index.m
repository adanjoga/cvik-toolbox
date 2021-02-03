function f = gd31index(X,clrs,K,pfun)
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
num = inf(K);
den = zeros(1,K);
for i = 1:K
    den(i) = max(max(DXX(clrs==i,clrs==i)));
    for j = setdiff(1:K,i)
        num(i,j) = mean2(DXX(clrs==i,clrs==j));
    end
end
f = min(num(:)/max(den));