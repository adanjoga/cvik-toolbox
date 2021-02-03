function D = cosdist(A,B)
% D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - mtimesx(A,'T',B);
D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - A'*B;