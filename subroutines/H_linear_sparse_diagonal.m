function H = H_linear(x)
% Evaluate the barrier for the nonnegative orthant
diag_vec = 1./(x.^2);
H = spdiags(diag_vec, 0, length(diag_vec), length(diag_vec));
end