function H = H_exp_tilde_dual_sparse_diagonal(input_vec)
% Compute the Hessian of the dual exp barrier (from the second pair)
% where N is the number of dual exp cones
N = length(input_vec)/3;
if floor(N)~=N
    display('Error: in g_exp_tilde_dual, the input vector has dimension not equal to 3*Ne');
end
H = sparse(3*N,3*N);
for k=1:N
    u = input_vec(3*k-2); v = input_vec(3*k-1); w = input_vec(3*k); r = log(-v/u); t = u-w+u*r;
    % Evaluate the Hessian of f_tilde_dual,
    % The (second) barrier function for the dual of the exponential cone
     H(3*k-2:3*k, 3*k-2:3*k) = [      1/(u*t) + r^2/t^2 + 1/u^2,             u*r/(v*t^2) - 1/(v*t),         -r/t^2;
                                (u*log(-v/u))/(v*t^2) - 1/(v*t), 1/v^2 + u/(v^2*t) + u^2/(v^2*t^2),     -u/(v*t^2);
                                                         -r/t^2,                        -u/(v*t^2),          1/t^2];
    
end

