function D = lapdist(A,B)
De = bsxfun(@plus,dot(B,B,1),dot(A,A,1)') - 2*(A'*B);%mtimesx(A,'T',B);
cosAB = 1 - 0.5*De;
D = (1-cosAB)./(1+cosAB);