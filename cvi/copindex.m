function f = copindex(X,clrs,K,pfun)
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
sumD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);          
end
intra = (sumD./Nk)';

% Dispersion inter-cluster
inter = zeros(1,K);
for i = 1:K
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XC = DXX(clrs==i,clrs==j);
        dd(:,k) = max(XC,[],2);
    end
    inter(i) = min(dd(:));
end
% Evaluacion del IVG
f = (1/N)*sum(Nk.*(intra./inter)');
