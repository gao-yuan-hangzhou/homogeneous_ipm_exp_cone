function grad = g_dual(z, dimension_info)
% Compute the gradient of the barrier function at z in K^*
% where K^* is a product of the nonnegative orthant, lorentz cones and dual exponential cones.
% x = [x_l; x_q; x_e], x_q = [x_q(1); ...; x_q(n_q)], x_e = [x_e(1); ...; x_e(n_e)];
if length(z)~=dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In g_dual(x, dimension_info), length of x does NOT match dimension_info!');
end
Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e;
% The linear part
g_l = g_linear(z(1:Nl));
% The (aggregated) quadratic part
g_q = zeros(sum(Nq),1); for k = 1:length(Nq) g_q(1+sum(Nq(1:k-1)):sum(Nq(1:k))) = g_lorentz(z(Nl+1+sum(Nq(1:k-1)):Nl+sum(Nq(1:k)))); end;
% The (aggregated) exponential part
g_e = g_exp_tilde_dual(z(end-3*Ne+1:end));
% concatenate all parts
grad = [g_l; g_q; g_e];
end

