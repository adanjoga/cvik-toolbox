function f = symdbindex(X,clrs,K,pfun)
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = Inf;
    return;
end
% Dispersion intra-cluster
M = NaN(K,size(X,2));
Si = zeros(1,K);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    
    dps = symdist(X(clrs == i,:),M(i,:),Nk(i),pfun);
    Si(i) = mean(dps);
end
% Dispersion inter-cluster
Mij = feval(pfun,M',M');% Distancia entre centroides
R = -inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        R(i,j) = (Si(i)+Si(j))/Mij(i,j);
    end
end
R = real(R); % Consider the real part only
% Sym DB index
f = mean(max(R,[],2));
