function Hessian = H_dual_NT(x, z, dimension_info)
% Hessian = blkdiag(H_linear_NT, H_soc_NT, H_exp_dual]
% where K^* is a product of the nonnegative orthant, lorentz cones and dual exponential cones.
% x = [x_l; x_q; x_e], x_q = [x_q(1); ...; x_q(n_q)], x_e = [x_e(1); ...; x_e(n_e)];

% Validate input
if length(z)~=dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In H_dual(z, dimension_info), length of x does NOT match dimension_info!');
end
% Obtain dimension information
Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e;
% Calculate H_l_NT
x_linear = x(1:Nl); z_linear = z(1:Nl);
diag_vec = z_linear ./ x_linear; H_l_NT = spdiags(diag_vec, 0, length(diag_vec), length(diag_vec));
% Calculate H_q_NT
H_q_NT = sparse(sum(Nq),sum(Nq)); for k = 1:length(Nq) H_q_NT(1+sum(Nq(1:k-1)):sum(Nq(1:k)), 1+sum(Nq(1:k-1)):sum(Nq(1:k))) = H_lorentz(z(Nl+1+sum(Nq(1:k-1)):Nl+sum(Nq(1:k)))); end;
% Calculate H_e
H_e = H_exp_tilde_dual_sparse_diagonal(z(end-3*Ne+1:end));
% Concatenate them diagonally
Hessian = blkdiag(H_l_NT, H_q_NT, H_e); 
end

