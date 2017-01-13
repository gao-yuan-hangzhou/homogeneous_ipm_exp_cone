function F = F_linear_sparse_diagonal(x)
% Find F such that F*F' = H_linear(x)
F = spdiags(1./x, 0, length(x), length(x));
end