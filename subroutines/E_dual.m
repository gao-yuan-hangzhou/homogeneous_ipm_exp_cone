function E = E_dual(z, dimension_info)
% Compute the Hessian of the barrier function at z in K^*,
% with the second-order cone part being replaced by its diagonal part
if length(z) ~= dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In H_dual(z, dimension_info), length of x does NOT match dimension_info!');
end

Nl = dimension_info.l; 
Nq = dimension_info.q; 
Ne = dimension_info.e;

% Calculate H_l, H_q and H_e
H_l = H_linear_sparse_diagonal(z(1:Nl));
H_q_diag = sparse(sum(Nq),sum(Nq)); 
for k = 1:length(Nq) 
    curr_indices = 1+sum(Nq(1:k-1)):sum(Nq(1:k));
    H_q_diag(curr_indices, curr_indices) = H_lorentz_diag_only(z(Nl+1+sum(Nq(1:k-1)):Nl+sum(Nq(1:k)))); 
end;
% disp(['density(H_q)=' num2str(nnz(H_q)/numel(H_q))]); 
H_e = H_exp_tilde_dual_sparse_diagonal(z(end-3*Ne+1:end));
E = sparse(blkdiag(H_l, H_q_diag, H_e)); % blkdiag(H_l, H_q, H_e) should already be sparse
end

