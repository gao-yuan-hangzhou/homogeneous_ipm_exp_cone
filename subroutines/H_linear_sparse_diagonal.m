function H = H_linear_sparse_diagonal(x)
% Evaluate the barrier for the nonnegative orthant
H = spdiags(1./(x.^2), 0, length(x), length(x));
end