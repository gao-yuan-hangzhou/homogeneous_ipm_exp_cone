function H_return = H_dual_NT(x, z, mu, dimension_info)
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
    w_q_curr = sqrt(gamma_Lorentz(z)/gamma_Lorentz(x));
    ksai_q_curr = z_q_curr/w_q_curr + w_q_curr*J_q_curr*x_q_curr;
    t_q_curr = sqrt(2)*J_q_curr*ksai_q_curr/(w_q_curr*gamma_Lorentz(ksai_q_curr));
    H_q(curr_q_idx, curr_q_idx) = -(1/(w_q_curr^2))*J_q_curr + t_q_curr*t_q_curr';
end;
H_e = H_exp_tilde_dual_sparse_diagonal(z(Nl+sum(Nq)+1:end));
H_return = blkdiag(H_l, H_q, mu*H_e); 
end

