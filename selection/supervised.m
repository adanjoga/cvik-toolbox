% Clustering Model Selection (CMS) using using a supervised metric (ari)
%   Clr is the best clustering solution from the PF in terms of ARI
function [Clr,ARIb,idx,ARIvalues] = supervised(CLRs, TLabels)

PFsize = size(CLRs,2);
ARIvalues = zeros(PFsize,1);
for i = 1:PFsize
    Yb = CLRs(:,i);
    ARIvalues(i) = pairwiseindex(TLabels,Yb);
end

[ARIb,idx] = max(ARIvalues,[],1);
Clr = CLRs(:,idx);
end