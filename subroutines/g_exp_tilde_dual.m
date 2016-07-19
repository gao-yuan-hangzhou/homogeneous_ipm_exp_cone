function grad = g_exp_tilde_dual(input_vec)
% Compute the gradient of the dual exponential barrier (second pair)
% N is the number of exp cones
Ne = length(input_vec)/3;
if floor(Ne)~=Ne
    display('Error: in g_exp_tilde_dual, the input vector z has dimension not equal to 3*Ne');
end
grad = zeros(3*Ne,1);
for k = 1:Ne
    u = input_vec(3*k-2); v = input_vec(3*k-1); w = input_vec(3*k); 
    r = log(-v/u); t = u-w+u*r;
    grad(3*k-2:3*k) = [-r/t-1/u; -1/v-u/(v*t); 1/t];
end
end

