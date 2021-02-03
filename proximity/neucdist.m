function D = neucdist(A, B)
A = bsxfun(@minus,A,mean(A,2));
B = bsxfun(@minus,B,mean(B,2));
A = bsxfun(@rdivide,A,std(A,0,2));
B = bsxfun(@rdivide,B,std(B,0,2));
% D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*mtimesx(A,'T',B));
D = sqrt(bsxfun(@plus,dot(B,B,1),dot(A,A,1)')- 2*(A'*B));