function f = silindex(X,clrs,K,pfun)
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
DX = DXX.^2;
S = zeros(1,K);
for i = 1:K
    XC = DX(clrs==clusts(i),clrs==clusts(i));
    ai = sum(XC,2)/max(Nk(i)-1, 1); % modificacion 25-10-16
    %ai = mean(XC,2);
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XK = DX(clrs==clusts(i),clrs==clusts(j));
        dd(:,k) = mean(XK,2);
    end
    bi = min(dd,[],2);
    S(i) = sum((bi-ai)./max([ai bi],[],2),1);
end
f = (1/N)*sum(S);
