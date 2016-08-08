addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);

% Paper reference: 
% https://github.com/gao-yuan-hangzhou/homogeneous_ipm_exp_cone/blob/master/test_Bernstein_approximation_chance_constrained_problem/note_PDF/bernstein_example.pdf

% Set the number of risky assets and number of underlying factors
n = 5; q = 2;

% Set alpha
alpha_risk = 0.05;

% Set the parameters
r0 = 1;
% eta(i) ~ LN(mu(i), sigma(i)^2), i = 1, ..., n_risky_assets
% zeta(l) ~ LN(v(l), theta(l)^2), l = 1, ..., n_factors
nv = zeros(q,1); theta = 0.1*ones(q,1);
gamma = abs(randn(n,q));
% Make sure 0<=rho<=0.1 and rho(1)<=...<=rho(n)
rho = sort(0.1*rand(n,1));
multiple = rho ./ (2*gamma*exp(nv+theta.^2/2));
for kk = 1:q
    gamma(:,kk) = gamma(:,kk) .* multiple;
end

log_E_eta = log(1 + rho/2); % E(eta(i)) = exp(mu(i)+sigma(i)^2/2), mu(i) = sigma(i)
mu = -1 + (2*log_E_eta+1).^(1/2); sigma = mu;

% Total number of random variables
d = n+q;

% Construct the discrete distributions
N = zeros(d,1);
for j = 1:n
    discrete_LN_vars{j} = xi_hat_discrete_LN(mu(j), sigma(j));
    N(j) = length(discrete_LN_vars{j}.vals);
end

for j = 1:q
    discrete_LN_vars{n+j} = xi_hat_discrete_LN(nv(j), theta(j));
    N(n+j) = length(discrete_LN_vars{n+j}.vals);
end

% Get the total number of discretized points
N_total_discretized_points = sum(N);

% Now the decision variables (including slack and auxiliary variables) are 

blk{1,1} = 'u'; blk{1,2} = 1;                                          % tau free
blk{2,1} = 'l'; blk{2,2} = n+2;                                        % x0, x(1), ... x(n), sx all >= 0               
blk{3,1} = 'u'; blk{3,2} = d+1;                                        % g0, g(1), ..., g(d) all free              
blk{4,1} = 'l'; blk{4,2} = 1;                                          % t0>=0
blk{5,1} = 'u'; blk{5,2} = d;                                          % s(1), ..., s(d) free
blk{6,1} = 'e'; blk{6,2} = 3*ones(N_total_discretized_points,1);       % [w(j,k); u(j,k); t(j,k)] in K_exp

% Total dimension of the decision vector
total_dim_dec_vec = 1+(n+2)+(d+1)+1+d+sum(blk{6,2});

% The objective function is 
% min -tau
c_cell{1} = -1;
for k = 2:6
    c_cell{k} = zeros(sum(blk{k,2}), 1);
end

% Get total number of constraints
m_total = (1+1+1) + d + N_total_discretized_points + d + (N_total_discretized_points - 1);
A = sparse(m_total, total_dim_dec_vec);
% The constraints are (code up the whole constraint matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 + x(1) + ... + x(n) + sx = 1 (1 constraint)
curr_row = 1;
A(curr_row,1+1:1+(2+n)) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g0 + s(1) + ... + s(n) - log_alpha*t + s0 = 0 (1 constraint)
curr_row = 2; 
idx_t = 1+(n+2)+(d+1)+1; idx_g0 = 1+(n+2)+1; idx_s0 = 1+(n+2)+(d+1)+2; idx_s_1_to_n = (1+(n+2)+(d+1)+2+1:1+(n+2)+(d+1)+2+1+(d-1));
A(curr_row,idx_g0) = 1; A(curr_row,idx_s_1_to_n) = 1; A(curr_row,idx_t) = -log(alpha_risk); A(curr_row,idx_s0) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g0 - tau + x0 = 0 (1 constraint)
curr_row = 3; 
idx_tau = 1; idx_x0 = 2;
A(curr_row,idx_g0) = 1; A(curr_row,idx_tau) = -1; A(3,idx_x0) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(j) + x(j) = 0, j = 1, ..., n (n constraints)
for j = 1:n
    curr_row = j+3;
    idx_gj = 1+(n+2)+1+j; idx_xj = 1+1+j;
    A(curr_row,idx_gj) = 1; A(curr_row,idx_xj) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(j) + gamma(:,j-n)*x(1:n) = 0, j = n+1, ..., n+q (q constraints)
for l = 1:q
    curr_row = n+3+l;
    idx_gl = 1+(n+2)+1+n+l; idx_x_1_to_n = (3:3+n-1);
    A(curr_row, idx_gl) = 1; A(curr_row,idx_x_1_to_n) = gamma(:,l)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w(j,k) - g(j)*v(j,k) + s(j) = 0, j = 1, ..., d, k = 1, ..., N(j) 
% (number of constraints = N_total_discretized_points)
for j = 1:d
    for k = 1:N(j)
        curr_row = (3+d) + sum(N(1:j-1)) + k;
        idx_wjk = (1+(n+2)+(d+1)+2+d) + 3*(sum(N(1:j-1))+k)-2; 
        idx_gj = 1+(n+2)+1+j; 
        idx_sj = 1+(n+2)+(d+1)+2+j;
        A(curr_row, idx_wjk) = 1;                                % w(j,k)
        A(curr_row, idx_gj) = -discrete_LN_vars{j}.vals(k);      % g(j)
        A(curr_row, idx_sj) = 1;                                 % s(j)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum(p(j,k)*u(j,k),k=1:N(j)) - t(j,1) = 0, any j (d constraints)
for j = 1:d
    curr_row = (3+d)+N_total_discretized_points+j;
    for k = 1:N(j)
        idx_ujk = 1+(n+2)+(d+1)+2+d+(3*(sum(N(1:j-1))+k)-1);
        A(curr_row, idx_ujk) = discrete_LN_vars{j}.prob_masses(k);
    end
    idx_tj1 = 1+(n+2)+(d+1)+2+d+3*(sum(N(1:j-1))+1);
    A(curr_row, idx_tj1) = -1;
end

% all t(j,k) are equal (number of constraints = N_total_discretized_points-1)
for j = 1:d
    for k = 1:N(j)
        if j==1 && k == 1 
            continue 
        end
        curr_row = (3+d)+N_total_discretized_points + d + (sum(N(1:j-1))+k-1);
        idx_t11 = 1+(n+2)+(d+1)+2+d+3;
        idx_tjk = 1+(n+2)+(d+1)+2+d+3*(sum(N(1:j-1))+k);
        A(curr_row, idx_t11) = 1;
        A(curr_row, idx_tjk) = -1;
    end
end

% Code up RHS vector b
b = [1; 0; 0; zeros(d,1); zeros(N_total_discretized_points,1); zeros(d,1); zeros(N_total_discretized_points-1,1)];
% 

% Split A into A_cell{kk}, kk = 1, ..., 6
A_cell{1} = A(:,blk{1,2});
A_cell{2} = A(:,blk{1,2}+1:blk{1,2}+blk{2,2});
A_cell{3} = A(:,blk{1,2}+blk{2,2}+1:blk{1,2}+blk{2,2}+blk{3,2});
A_cell{4} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2});
A_cell{5} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2});
A_cell{6} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2}+sum(blk{6,2}));
[obj_val, x_re, y_re, z_re, info] = hsd_lqeu(blk, A_cell, c_cell, b);