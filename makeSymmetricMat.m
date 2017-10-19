function A = makeSymmetricMat(A)
At = A';
idx1 = (A~=At);
idx2 = (A==0);
idx = (idx1&idx2);
idx3 = find(idx>0);
A(idx3) = At(idx3);
end