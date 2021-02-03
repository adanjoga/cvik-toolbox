function f = sfindex(X,clrs,K,pfun)
global Xmean;

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
wcd = sum(sumD./Nk);

% Dispersion inter-cluster
DMMg = feval(pfun,M',Xmean);
bcd = sum(Nk.*DMMg)/(N*K);

% Evaluacion del IVG
f = 1 - (1/(exp(exp(bcd-wcd))));