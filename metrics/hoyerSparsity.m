function S = hoyerSparsity(w)
%HOYERSPARSITY Hoyer sparsity measure normalized to [0,1].
%
% S = 0 -> dense/equal-magnitude vector
% S = 1 -> maximally sparse vector (one nonzero coefficient)

M = numel(w);

l1 = sum(abs(w));
l2 = norm(w,2);

if l2 == 0
    S = 1;
else
    S = (sqrt(M) - l1/l2) / (sqrt(M)-1);
    S = min(max(S,0),1);
end
end
