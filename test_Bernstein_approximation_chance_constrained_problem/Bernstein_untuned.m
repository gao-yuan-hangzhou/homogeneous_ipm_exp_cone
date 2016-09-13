addpath(fileparts(pwd)); 
addpath([fileparts(pwd), '/subroutines']);
addpath([fileparts(pwd), '/subroutines']);
clear;
% Whether to use CVX to solve the Bernstein approximation
is_using_cvx = false;

% Paper reference: 
% https://github.com/gao-yuan-hangzhou/homogeneous_ipm_exp_cone/blob/master/test_Bernstein_approximation_chance_constrained_problem/note_PDF/bernstein_example.pdf

% Set the number of risky assets and number of underlying factors
n = 64; 
q = 8; 
disp(['number_of_risky_assets = ' num2str(n)]);
disp(['number_of_common_factors = ' num2str(q)]);

% Set alpha
alpha_risk = 0.05;

% Set the parameters
r0 = 1;
% eta(i) ~ LN(mu(i), sigma(i)^2), i = 1, ..., n_risky_assets
% zeta(l) ~ LN(v(l), theta(l)^2), l = 1, ..., n_factors
nu = 0.01*randn(q,1); theta = rand(q,1);
gamma = abs(randn(n,q));
% Make sure 0<=rho<=0.1 and rho(1)<=...<=rho(n)
rho = sort(0.1*rand(n,1));
multiple = (rho/2) ./ (gamma*exp(nu+theta.^2/2));
for kk = 1:q
    gamma(:,kk) = gamma(:,kk) .* multiple;
end

log_E_eta = log(1 + rho/2); % E(eta(i)) = exp(mu(i)+sigma(i)^2/2), mu(i) = sigma(i)
mu = -1 + (2*log_E_eta+1).^(1/2); sig = mu; mu = log_E_eta - sig.^2/2;

% Save the parameters for debugging
% save('n', 'q', 'nu', 'theta', 'mu', 'sig', 'rho', 'gamma');
% load('n', 'q', 'nu', 'theta', 'mu', 'sig', 'rho', 'gamma');
%load('text_example_CVX_failed_2.mat', 'n', 'q', 'nu', 'theta', 'mu', 'sig', 'rho', 'gamma');

% Total number of random variables
d = n+q;

% Construct the discrete distributions
eps_th = 1e-4;
Del_resol = 0.01;
N = zeros(d,1);
for j = 1:n
    discrete_LN_vars{j} = xi_hat_discrete_LN(mu(j), sig(j), eps_th, Del_resol);
    v{j} = discrete_LN_vars{j}.vals; 
    p{j} = discrete_LN_vars{j}.prob_masses;
    N(j) = length(v{j});
end

for j = 1:q
    discrete_LN_vars{n+j} = xi_hat_discrete_LN(nu(j), theta(j), eps_th, Del_resol);
    v{n+j} = discrete_LN_vars{n+j}.vals; 
    p{n+j} = discrete_LN_vars{n+j}.prob_masses;
    N(n+j) = length(v{n+j});
end

t_begin_hsd_lqeu = cputime;
% Get the total number of discretized points
N_total = sum(N);

% Now the decision variables (including slack and auxiliary variables) are 
blk{1,1} = 'u'; blk{1,2} = 1;                                          % tau free
blk{2,1} = 'l'; blk{2,2} = n+2;                                        % x0, x(1), ... x(n), sx all >= 0               
blk{3,1} = 'u'; blk{3,2} = d+1;                                        % g0, g(1), ..., g(d) all free              
blk{4,1} = 'l'; blk{4,2} = 1;                                          % t0>=0
blk{5,1} = 'u'; blk{5,2} = d;                                          % s(1), ..., s(d) free
blk{6,1} = 'e'; blk{6,2} = 3*ones(N_total,1);                          % [w(j,k); u(j,k); t(j,k)] in K_exp

% Total dimension of the decision vector
total_dim = 1+(n+2)+(d+1)+1+d+3*N_total;

% The objective function is 
% min -tau
c_cell{1} = -1;
for k = 2:6
    c_cell{k} = zeros(sum(blk{k,2}),1);
end

% Get total number of constraints
m_total = (1+1+1) + (n+q) + N_total + d + N_total;
A = sparse(m_total, total_dim);

% Set some indices of the decision variables
idx_tau = 1; idx_x0 = 2; idx_x1 = 3; idx_x1_to_xn = (3:n+2); 
idx_x0_to_sx = (2:n+3); idx_sx = n+3; idx_g0 = n+4; idx_g_1 = n+5;
idx_g_1_to_d = (n+5:n+5+(d-1)); idx_t0 = n+5+(d-1)+1; 
idx_s1 = n+5+(d-1)+2; idx_s1_to_sn = n+5+(d-1)+2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The constraints are (code up the whole constraint matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 + x(1) + ... + x(n) + sx = 1 (1 constraint)
curr_row = 1;
A(curr_row,idx_x0_to_sx) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g0 + (s(1) + ... + s(d)) - log_alpha*t0 = 0 (1 constraint)
curr_row = 2; 
idx_s_1_to_d = (1+(n+2)+(d+1)+1+1:1+(n+2)+(d+1)+1+d);
A(curr_row,idx_g0) = 1;
A(curr_row,idx_s_1_to_d) = 1; 
A(curr_row,idx_t0) = -log(alpha_risk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g0 - tau + x0 = 0 (1 constraint)
curr_row = 3; 
idx_tau = 1; 
idx_x0 = 2;
A(curr_row,idx_g0) = 1; 
A(curr_row,idx_tau) = -1; 
A(curr_row,idx_x0) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(j) + x(j) = 0, j = 1, ..., n (n constraints)
for j = 1:n
    curr_row = 3+j;
    idx_gj = 1+(n+2)+1+j; 
    idx_xj = 1+1+j;
    A(curr_row,idx_gj) = 1; 
    A(curr_row,idx_xj) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g(n+l) + gamma(:,l)*x(1:n) = 0, l = 1, ..., q (q constraints)
for l = 1:q
    curr_row = 3+n+l;
    idx_gj = 1+(n+2)+1+n+l; 
    idx_x1_to_xn = (3:2+n);
    A(curr_row,idx_gj) = 1; 
    A(curr_row,idx_x1_to_xn) = gamma(:,l)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w(j,k) - g(j)*v(j,k) + s(j) = 0, j = 1, ..., d, k = 1, ..., N(j) 
% (number of constraints = N_total_discretized_points)
A6 = sparse(N_total, total_dim);
for j = 1:d
    for k = 1:N(j)
        curr_row = 3 + d + (sum(N(1:j-1))+k);
        idx_wjk = (1+(n+2)+(d+1)+1+d) + (3*(sum(N(1:j-1))+k)-2);
        idx_gj = 1+(n+2)+1+j;
        idx_sj = 1+(n+2)+(d+1)+1+j;
        A(curr_row,idx_wjk) = 1;                               % w(j,k)
        A(curr_row,idx_gj) = -v{j}(k);                         % g(j)
        A(curr_row,idx_sj) = 1;                                % s(j)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sum(p(j,k)*u(j,k),k=1:N(j)) - t(j,1) = 0, any j (d constraints)
A7 = sparse(d, total_dim);
for j = 1:d
    curr_row = 3 + d + N_total+j;
    for k = 1:N(j)
        idx_ujk = 1+(n+2)+(d+1)+1+d+(3*(sum(N(1:j-1))+k)-1);
        A(curr_row,idx_ujk) = p{j}(k);
        A7(j,idx_ujk) = p{j}(k);
    end
    idx_t0 = 1+(n+2)+(d+1)+1;
    A(curr_row,idx_t0) = -1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t0 - t(j,k) = 0, j = 1, ..., d, k = 1, ..., N(j)
for j = 1:d
    for k = 1:N(j)
        curr_row = 3 + d + N_total + d + (sum(N(1:j-1))+k);
        idx_t0 = 1+(n+2)+(d+1)+1;
        idx_tjk = 1+(n+2)+(d+1)+1+d+3*(sum(N(1:j-1))+k);
        A(curr_row,idx_t0) = 1;
        A(curr_row,idx_tjk) = -1;
    end
end

% Code up RHS vector b
b = [1; 0; 0; zeros(d,1); zeros(N_total,1); zeros(d,1); zeros(N_total,1)];

% Split A into A_cell{kk}, kk = 1, ..., 6
A_cell{1} = A(:,blk{1,2});
A_cell{2} = A(:,blk{1,2}+1:blk{1,2}+blk{2,2});
A_cell{3} = A(:,blk{1,2}+blk{2,2}+1:blk{1,2}+blk{2,2}+blk{3,2});
A_cell{4} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2});
A_cell{5} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2});
A_cell{6} = A(:,blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2}+1:blk{1,2}+blk{2,2}+blk{3,2}+blk{4,2}+blk{5,2}+sum(blk{6,2}));

% Display the time for constructing the input
disp('Done constructing blk, A_cell, c_cell, b for hsd_sqeu_Schur!');

% ============== Call the solver ==============
[obj_val, x_re, y_re, z_re, info] = hsd_lqeu_fast_lu(blk, A_cell, c_cell, b);
[obj_val, x_re, y_re, z_re, info] = hsd_lqeu(blk, A_cell, c_cell, b);
% =============================================
obj_tau_minus_one = -obj_val(2) - 1;
hsd_lqeu_x_opt = x_re{2}(2:end-1);
t_end_hsd_lqeu = cputime;
t_elapsed_hsd_lqeu = t_end_hsd_lqeu - t_begin_hsd_lqeu;

if is_using_cvx
    % Solve the Bernstein approximation
    cvx_clear
    t_begin_cvx = cputime;
    cvx_begin
    log_alpha = log(alpha_risk);
    variables tau x0 x(n) g0 g(d) s(d) t w(N_total) u(N_total)

    maximize(tau-1)
        x0 >= 0; x>=0; 
        x0 + sum(x) <= 1;
        g0 + sum(s) - log_alpha*t == 0;
        g0 == tau - x0;
        for j = 1:n 
            g(j) == -x(j); 
        end;
        for l = 1:q
            g(n+l) == -gamma(:,l)'*x;
        end
        for j = 1:d
            % sum(p(j,k) * u(j,k), k = 1:N(j)) = t
            curr_idx_u = (sum(N(1:j-1))+1:sum(N(1:j)));
            p{j}'*u(curr_idx_u) == t;
            for k = 1:N(j)
                curr_idx_jk = sum(N(1:j-1)) + k;
                w(curr_idx_jk) == v{j}(k)*g(j) - s(j);
                % Note that CVX defines exponential cone in a slight different way
                % by interchanging the role of y and z in (x,y,z) in K_exp
                {w(curr_idx_jk), t, u(curr_idx_jk)} <In> exponential;
            end
        end
    cvx_end
    t_end_cvx = cputime;
    t_elapsed_cvx = t_end_cvx - t_begin_cvx;
    obj_tau_minus_one_cvx = cvx_optval;
    cvx_x_opt = cvx_value(x);
end
    
% Solve for the nominal dterministic optimal value (all random variables replaced by their means)
cvx_clear
cvx_begin
variables tau x0 x(n)
maximize(tau-1)
    tau <= (1+rho)'*x;
    sum(x)<=1; x>=0;
cvx_end

opt_nominal = cvx_optval;
disp(['number_of_risky_assets = ' num2str(n)]);
disp(['number_of_common_factors = ' num2str(q)]);
disp(['risk tolerance alpha = ' num2str(alpha_risk)]); disp(' ');
disp(['nominal optimal value = ' num2str(opt_nominal)]);
disp(['maximum return by hsd_lqeu_Schur = ' num2str(obj_tau_minus_one)]);
disp(['running time = ' num2str(t_elapsed_hsd_lqeu) ' seconds']);
disp(['optimal risky asset allocation vector obtained by hsd_lqeu_Schur (largest 10 entries) = ' num2str(obj_tau_minus_one)]); 
vec1 = sort(hsd_lqeu_x_opt, 'descend'); disp(vec1(1:min(10,end))');
if is_using_cvx
    disp(['maximum return by CVX (SOCP approximation and calling SDPT3) = ' num2str(obj_tau_minus_one_cvx)]);
    disp(['running time = ' num2str(t_elapsed_cvx) ' seconds']);
    disp('optimal risky asset allocation vector obtained by CVX (largest 10 entries) = '); 
    vec2 = sort(cvx_x_opt, 'descend'); disp(vec2(1:min(10,end))');
end