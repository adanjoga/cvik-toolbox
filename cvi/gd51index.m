function f = gd51index(X,clrs,K,pfun)
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
M = NaN(K,size(X,2));
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);  
end

num = inf(K);
den = zeros(1,K);
for i = 1:K
    den(i) = max(max(DXX(clrs==i,clrs==i)));
    for j = setdiff(1:K,i)
        DxS = feval(pfun,X(clrs==i,:)',M(j,:)');
        DyT = feval(pfun,X(clrs==j,:)',M(i,:)');
        num(i,j) = (1/(Nk(i)+Nk(j)))*(sum(DxS)+sum(DyT));
    end
end
f = min(num(:)/max(den));