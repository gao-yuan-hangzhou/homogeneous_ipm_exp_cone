function output = is_in_product_cone(x_bar, dimension_info)
% Check whether x_bar = [x; y; z; tau; kappa; theta] satisfies
% x in int(K), z in int(K^*), kappa > 0, tau > 0
output = true;
% Find the dimensions of the linear, quadratic and exponential parts
Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e; Nt = Nl + sum(Nq) + 3*Ne; m= dimension_info.m;
% Extract x, y, z, tau, kappa
x = x_bar(1:Nt); z = x_bar(Nt+m+1:2*Nt+m); % y = x_bar(Nt+1:Nt+m);
tau = x_bar(2*Nt+m+1); kappa = x_bar(2*Nt+m+2); % theta = x_bar(2*Nt+m+3);
% Further extract x_l, x_q, x_e, z_l, z_q, z_e_dual
x_l = x(1:Nl); x_q = x(Nl+1:end-Ne); x_e = x(end-3*Ne+1:end);
z_l = z(1:Nl); z_q = z(Nl+1:end-Ne); z_e_dual = z(end-3*Ne+1:end);

% Check tau, kappa, x_l, z_l, x_e, z_e_dual
if ~(tau>0) || ~(kappa>0) || ~(prod(x_l>0)) || ~(prod(z_l>0)) || ~is_in_exp_cone_interior(x_e) || ~is_in_dual_exp_cone_interior(z_e_dual)
    output = false; 
end
% Check x_q, z_q
for k = 1:length(Nq)
    if ~(is_in_lorentz_cone_interior(x_q(sum(Nq(1:k-1))+1:sum(Nq(1:k))))) || ~(is_in_lorentz_cone_interior(z_q(sum(Nq(1:k-1))+1:sum(Nq(1:k)))))
        output = false; break;
    end
end
end

