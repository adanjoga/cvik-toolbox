% Point Symmetry-Distance
function dps = symdist(Xk,Mk,Nk,pfun)
m2 = 2*Mk;
m2x = m2(ones(1,Nk),:)-Xk;
de = feval(pfun,m2x',Xk');
de = sort(de,2);
dps = mean(de(:,1:2),2);