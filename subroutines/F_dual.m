function F = F_dual(z, dimension_info)
% F is a nonsingular matrix such that F * F' = H_dual
% Compute the Hessian of the barrier function at z in K^*,
% where K^* is a product of the nonnegative orthant, lorentz cones and dual exponential cones.
% x = [x_l; x_q; x_e], x_q = [x_q(1); ...; x_q(n_q)], x_e = [x_e(1); ...; x_e(n_e)];
if length(z)~=dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In H_dual(z, dimension_info), length of x does NOT match dimension_info!');
end
Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e;

% Calculate H_l, H_q and H_e
F_l = F_linear_sparse_diagonal(z(1:Nl));
F_q = sparse(sum(Nq),sum(Nq)); 
for k = 1:length(Nq) 
    curr_indices = 1+sum(Nq(1:k-1)):sum(Nq(1:k));
    F_q(curr_indices, curr_indices) = F_lorentz(z(Nl+1+sum(Nq(1:k-1)):Nl+sum(Nq(1:k)))); 
end;
% disp(['density(H_q)=' num2str(nnz(H_q)/numel(H_q))]); 
F_e = F_exp_tilde_dual_sparse_diagonal(z(end-3*Ne+1:end));
F = sparse(blkdiag(F_l, F_q, F_e)); % blkdiag(H_l, H_q, H_e) should already be sparse
end
