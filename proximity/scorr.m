function D = scorr(~,B)
global Xrank;
A = Xrank; % A = ranking(A);
B = ranking(B);
A = bsxfun(@minus,A,mean(A,2));
B = bsxfun(@minus,B,mean(B,2));
% D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - mtimesx(A,'T',B);
D = bsxfun(@plus,0.5*dot(B,B,1),0.5*dot(A,A,1)') - A'*B;