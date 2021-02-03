function Vi = boundConstraint(Vi,K,d,Xmin,Xmax)
A = Xmin(ones(1,K),:);
B = Xmax(ones(1,K),:);
R = A+(B-A).*rand(K,d);
idx = (Vi>B)|(Vi<A);
Vi(idx) = R(idx);