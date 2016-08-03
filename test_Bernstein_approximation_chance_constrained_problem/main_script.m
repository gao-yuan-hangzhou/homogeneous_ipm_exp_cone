% Bernstein approximation vs Scenario approach on a toy example
% as in the paper of Nemirovski and Shapiro

% Set the number of risky assets and number of underlying factors
n = 64; q = 4;

% Set the parameters
r0 = 1;
% eta(i) ~ LN(mu(i), sigma(i)^2), i = 1, ..., n_risky_assets
% zeta(l) ~ LN(v(l), theta(l)^2), l = 1, ..., n_factors
v = zeros(q,1); theta = 0.1*ones(q,1);
gamma = abs(randn(n, q));
rho = 2*gamma * exp(v+theta.^2/2);

log_E_eta = log(1 + rho/2); % E(eta(i)) = exp(mu(i)+sigma(i)^2/2), mu(i) = sigma(i)
mu = -1 + (2*log_E_eta+1).^(1/2); sigma = mu;

% The original problem is
% max{(tau - 1)|x>=0, tau free} 
% s.t. P(tau>sum(r(j)x(j),j=0:n_risky_assets))<=alpha, sum(x(j),j=0:n_risky_assets)<=1
% The problem is equivalent to 
% max{(tau - 1)|x,tau>=0} 
% s.t. P(F(xBar,xi)<=0)>=1-alpha, where xBar = [tau; x0; x(1:n)], xi = [eta(1:n); zeta(1:q)]
% F(xBar,xi) = g_0(x) + sum(eta(j)*xg_j(x),j=1:n) + sum(zeta(l)*h_l(x),l=1:q), 
% where g_j(xBar) = -x(j), h_l(xBar) = -sum(gamma(j,l)*x(j),j=1:n)

% The Bernstein approximation of the above problem is (converted into a minimization problem)
% min{-tau+1|xBar>=0}
% s.t. g_0(xBar)+sum(t*logMGF_j(t^(-1)*g_j(xBar)),j=1:d)+s <= 0, t>=0, s >= 0, x>=0

% xBar = [tau; x0; x(1:n); s0; w(k,j); ]
