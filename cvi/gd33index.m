function f = gd33index(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -inf;
    return;
end

% Dispersion intra-cluster
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)'));          
end
den = 2*(sumD./Nk);

% Dispersion inter-cluster
num = inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        num(i,j) = mean2(DXX(clrs==i,clrs==j));
    end
end
% Evaluacion del IVG
f = min(num(:)/max(den));