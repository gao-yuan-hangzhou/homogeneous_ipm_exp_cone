function discrete_pmf = xi_hat(mu, sigma, eps_th, Del_resol)
% For xi ~ LN(mu, sigma^2), generate xi_hat represented by [v(1:N), p(1:N)] where 
% p(i) = P(xi_hat = v(i)) and p(1)+...+p(N) = 1, which is a 
% discrete approximation of the lognormal random variable xi.
% The apprixmation scheme ensures that 0<xi_hat <= xi, and 
% P(xi_hat/xi >= exp(-Del_resol)) >= 1 - eps_th.

% Return a structure array containing vals and prob_masses
% such that P(xi = vals(k)) = prob_masses(k)

% Find R such that  normcdf(+inf) - normcdf(R) = threshold_prob/2
% and construct the "bins" and assign values of b(k)
R = norminv(1-eps_th/2);
bins_length = Del_resol/sigma;
a = [(-R:bins_length:R)'; R];
b0 = 0; b = exp(sigma*a+mu);
v_mass0 = normcdf(a(1)); 
v_mass = zeros(length(a),1);
for k = 1:length(v_mass)-1
    v_mass(k) = normcdf(a(k+1)) - normcdf(a(k));
end
v_mass(end) = 1 - normcdf(a(end));

discrete_pmf.vals = [b0; b];
discrete_pmf.prob_masses = [v_mass0; v_mass];

end

