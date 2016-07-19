function H = H_linear(x)
% Evaluate the barrier for the nonnegative orthant
H = diag(1./(x.^2));
end