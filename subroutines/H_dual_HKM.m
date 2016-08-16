function H_return = H_dual_HKM(x, z, mu_xz, dimension_info)
% Compute the Hessian of the barrier function at w (NT scaling point) in the 
% relative interior of the product cone K^*

if length(z)~=dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In H_dual(z, dimension_info), length of x does NOT match dimension_info!');
end

Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e;

% Calculate H_l(wl)
H_l = H_linear_sparse_diagonal(x(1:Nl)./z(1:Nl));
H_q = sparse(sum(Nq),sum(Nq)); 
for k = 1:length(Nq)
    Nq_curr = Nq(k);
    J_q_curr = diag([1;-1*ones(Nq_curr-1,1)]);
    curr_q_idx = (1+sum(Nq(1:k-1)):sum(Nq(1:k)));
    x_q_curr = x(Nl+curr_q_idx);
    z_q_curr = z(Nl+curr_q_idx);
    term1 = - x_q_curr'*z_q_curr/(z_q_curr(1)^2 - z_q_curr(2:end)'*z_q_curr(2:end))*J_q_curr;
    term2 = x_q_curr * Lorentz_inv(z_q_curr)';
    term3 = Lorentz_inv(z_q_curr) * z_q_curr';
    H_q(curr_q_idx, curr_q_idx) = term1 + term2 + term3;
end;
% disp(['density(H_q)=' num2str(nnz(H_q)/numel(H_q))]); 
H_e = H_exp_tilde_dual_sparse_diagonal(z(end-3*Ne+1:end));
H_return = blkdiag(H_l, H_q, mu_xz*H_e); 
end

