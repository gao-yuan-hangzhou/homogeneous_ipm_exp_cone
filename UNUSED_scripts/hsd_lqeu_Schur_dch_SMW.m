function [obj_val, x_return,y_return,z_return, result_info] = hsd_lqeu_Schur_dense_column(blk, A_cell, c_cell, b, rel_eps, max_iter_count)
format long;
addpath([fileparts(pwd), '/subroutines']);
disp('=== hsd_lqeu_Schur_dch_SMW started... ===');
disp('This program uses Sherman-Morrison-Woodbury formula in solving the Schur complement equation...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves problems of the following form
% min sum_j (c(j)'x(j)) s.t. sum_j (A(j)x(j)) = b, x(j) in K(j) (or free), j=1,2,...,N
% where K(j) can be a nonnegative orthant of arbitrary dimension, a Lorentz cone, product of Lorentz
% cones, product of exponential cones or just R^n. See below for the definitions.
% Definitions:
% R_+ = {x: x>=0}.
% Q(n) = {x in R^n: x(1)>=||x(2:n)||}, where || || is the Euclidean 2-norm.
% K_exp = closure{(x1,x2,x3): x2>=0, x3>0, exp(x1/x3)<=x2/x3} is the exponential cone.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The exact input format is described as follows.
% [obj,X,y,Z,info] = hsd_lqe(blk,At,c,b,rel_eps)
% Input arguments:
% blk: a cell array describing the dimensions of input data
% A_cell: the data for A in Ax=b
% c_cell: the data for c in min c'x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: Suppose one has 
% x(1) in (R_+)^15, x(2) in Q(3)�Q(4)�Q(6), x(3) in (K_exp)^2, x(4) in Q(8), 
% x(5) in (K_exp)^7, x(6) in R^9, m=23.
% Then necessarily c = (c{1}, c{2}, c{3}, c{4}, c{5}) is a vector of dimension 
% 15 + (3+4+6) + 3*2 + 8 + 3*7 + 9 = 72
% b is a vector of dimension m=23.
% The input arguments are
% blk{1,1} = 'l'; blk{1,2} = 15;                      At{1} is a sparse 23-by-15 matrix
% blk{2,1} = 'q'; blk{2,2} = [3; 4; 6];               At{2} is a sparse 23-by-13 matrix  (3+4+6=13)
% blk{3,1} = 'e'; blk{3,2} = [3;3];                   At{3} is a sparse matrix of dimension 23-by-6
% blk{4,1} = 'q'; blk{4,2} = [8];                     At{4} is a sparse matrix of dimension 23-by-8
% blk{5,1} = 'e'; blk{5,2} = [3;3;3;3;3;3;3];         At{5} is a sparse matrix of dimension 23-by-21
% blk{6,1} = 'u'; blk{6,2} = 9;                       At{6} is a sparse 23-by-9 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All subroutines are in the folder "subroutines"
% Get the current CPU time at the beginning
t_begin = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First of all, convert all 'u' into 'l' using the simplest scheme x = x_pos - x_neg
% and keep track of the original indices of 'u'
indices_of_u = zeros(size(blk,1),1);
for k = 1:size(blk,1)
    if blk{k,1} == 'u'
        indices_of_u(k) = 1;
        blk{k,1} = 'l'; blk{k,2} = 2*blk{k,2};
        A_cell{k} = [A_cell{k}, -A_cell{k}]; 
        c_cell{k} = [c_cell{k}; -c_cell{k}];
    end
end

if sum(indices_of_u) > 0
    disp('Unrestricted block(s) detected and converted into linear block(s)!');
end

% Obtain Nl, Nq = [q(1), ..., q(n_1)], Ne (number of exponential cones) and m
Nl = 0; Nq = []; Ne = 0; m = length(b);
for k = 1:size(blk,1)
    if blk{k,1} == 'l' Nl = Nl+blk{k,2}; 
    elseif blk{k,1} == 'q' Nq = [Nq; blk{k,2}];
    elseif blk{k,1} == 'e' Ne = Ne + length(blk{k,2});
    else display('Error: blk{k,1} is not one of l, q or e!');
    end;
end;
% Compute the total dimension of x = [xl; xq; xe]
dim_x = Nl + sum(Nq) + 3*Ne;
% Define v, the parameter for the logarithmic barrier for the product cone K
v_param = Nl + length(Nq) + 3*Ne;
% Define structure dimension_info that contains Nl, Nq and Ne
dimension_info.l = Nl; dimension_info.q = Nq; dimension_info.e = Ne; dimension_info.m = m;
% Construct A = [Al, Aq, Ae]
Al = sparse(m,0); Aq = sparse(m,0); Ae = sparse(m,0);
for k=1:size(blk,1)
    if blk{k,1} == 'l' 
        Al = [Al, A_cell{k}];
    elseif blk{k,1} == 'q' 
        Aq = [Aq, A_cell{k}];
    elseif blk{k,1} == 'e' 
        Ae = [Ae, A_cell{k}];
    else
        disp(['Error: blk{' num2str(k) ',1} is not one of l, q or e!']);
    end;
end;
A = [Al, Aq, Ae];

% Contruct c = [cl; cq; ce]
cl = sparse(0,1); cq = sparse(0,1); ce = sparse(0,1);
for k=1:size(blk,1)
    if blk{k,1} == 'l' 
        cl = [cl; c_cell{k}];
    elseif blk{k,1} == 'q' 
        cq = [cq; c_cell{k}];
    elseif blk{k,1} == 'e' 
        ce = [ce; c_cell{k}];
    else display('Error: blk{k,1} is not one of l, q or e!');
    end;
end;
c = [cl; cq; ce];

% Check whether Ax=b has a solution
% disp('Check whether A has full row rank and Ax=b has a solution...');
% [~, U_A] = lu(A); [~,U_full] = lu([A,b]);
% if sprank(U_A) < sprank(U_full)
%     display('Error: Ax=b has no solution! The primal problem is infeasible.'); return;
% elseif sprank(U_A) < m
%     display('Error: A is not FULL ROW RANK!');
%     return;
% else
%     disp('Ok: A has full row rank and Ax = b');
% end

% Now initialize (x, y, z, tau, kappa, theta)
exp_cone_center = [-1.0151; 1.2590; 0.5560];
id_l = ones(Nl,1);
id_q = zeros(sum(Nq),1); 
for k = 1:length(Nq) 
    id_q(sum(Nq(1:k-1))+1) = 1; 
end;
id_e = repmat(exp_cone_center,Ne,1);
x0 = [id_l; id_q; id_e]; x = x0; y0 = zeros(m,1); y = y0; z0 = x0; z = z0;
tau = 1; kappa = 1; theta0 = 1.0; theta = theta0;
% Compute the auxiliary parameters which completely define (HSD)
b_bar = (b*tau - A*x)/theta; c_bar = (c*tau - A'*y - z)/theta;
g_bar = (c'*x - b'*y + kappa)/theta; alpha_bar = (x'*z + tau*kappa)/theta;

% Set the default relative accuracy
if nargin < 5
    rel_eps = 1e-8;
end

result_info.solution_status = 'undetermined';

% Set the default max number of iterations
if nargin < 6
    max_iter_count = 500;
end

% The main loop
for it_count =1:max_iter_count
    % disp(['(<x,z> + tau*kappa) - alpha_bar*theta = ' num2str(x'*z+tau*kappa - alpha_bar*theta)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the HKM-like search direction pred_dir = [dx; dy; dz; dtau; dkappa]
    % as in SDPT3 Guide: http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf
    A_hat = [A; -c'; c_bar']; % A_hat is (m+2) by dim_x
    y_hat = [y; tau; theta]; % y_hat has dimension (m+2)
    B_hat = [sparse(m,m), -b, b_bar; b', 0, g_bar; -b_bar', -g_bar, 0]; % B_hat is (m+2) by (m+2)
    % mu_xz = x'*z/v; mu_hat = theta;
    mu_hat = theta; % mu_hat = (x'*z +tau*kappa)/(v_param+1);
    Rp_hat = [sparse(m,1); kappa; -alpha_bar] - A_hat*x - B_hat*y_hat;
    Rd = -A_hat'*y_hat - z;
    Rc = -x; % sigma = 0 for predictor direction system, as explained below
    Rt = - kappa; % Rt = mu_hat/tau - kappa;
    H = H_dual(z, dimension_info); H = mu_hat*H; % to fix Matlab's memory allocation issue
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve equation (28) and then (27) with sigma = 0 for the predictor direction
    % Make use of the dense column handling techniques described in SDPT3 guide
    % LHS_Schur_comp_eq = A_hat*H*A_hat' + B_hat + diag([zeros(m,1); kappa/tau; 0]);
    h_hat = sparse(Rp_hat + A_hat*(H*Rd-Rc) + [sparse(m,1); Rt; 0]);
    Msp = A_hat*H*A_hat'+diag([sparse(m,1); kappa/tau; 0]);
    % disp(['nnz(R) = ' num2str(nnz(R)) ', density(R) = ' num2str(nnz(R)/numel(R))]);    
    % disp(['nnz(Msp) = ' num2str(nnz(Msp)) ', density(Msp) = ' num2str(nnz(Msp)/numel(Msp))]);
    U1 = [-b; 0; -g_bar/2]; U2 = [b_bar; g_bar/2; 0];
    V1 = [zeros(m,1); 1; 0]; V2 = [zeros(m+1,1); 1];
    U = [U1, U2, V1, V2];
    % D = [zeros(2), eye(2); -eye(2), zeros(2)];
    D_inv = [zeros(2), -eye(2); eye(2), zeros(2)];
    % Now solve for dy_hat using Sherman-Morrison formula
    % Calculate inv(Msp)*h_hat = R\(R'\h_hat)
    Msp_inv_h_hat = Msp\h_hat;
    % Calculate the 4-by-4 matrix G = D_inv + U'*inv(Msp)*U
    Msp_inv_U = Msp\U;
    G = D_inv + U'*Msp_inv_U;
    % Implicity calculate P = eye(m+2) - inv(Msp)*U*G_inv*U' and then y_hat = P*w_hat
    % P = speye(m+2) - Msp_inv_U*(G\U');
    % Calculate dy_hat_pred = LHS_Schur_comp_eq \ h_hat; % (28)
    dy_hat_pred = Msp_inv_h_hat - Msp_inv_U*(G\(U'*Msp_inv_h_hat));
    % Solve for dy_pred using the naive M\h_hat for DEBUGGING
    % M_schur = Msp+U*D*U';
    % dy_hat_pred = M_schur\h_hat;
    % disp(['error in dy_hat_pred = ' num2str(norm(dy_hat_pred-dy_hat_pred_sm)/norm(dy_hat_pred))]);
    % dy_hat_pred = dy_hat_pred_sm;
    dy_pred = dy_hat_pred(1:m);
    dtau_pred = dy_hat_pred(m+1);
    dtheta_pred = dy_hat_pred(m+2);
    dz_pred = Rd - A_hat'*dy_hat_pred; % 2nd equation of (27)
    dx_pred = Rc - H * dz_pred; % 3rd equation of (27)
    dkappa_pred = Rt - (kappa/tau)*dtau_pred; % 4th equation of (27)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute alpha_p (approximately)
    x_bar = [x; y; z; tau; kappa; theta];
    pred_dir = [dx_pred; dy_pred; dz_pred; dtau_pred; dkappa_pred; dtheta_pred];
    alpha_p = find_alpha_max(x_bar, pred_dir, dimension_info);
    % Set sigma (parameter for the system for the combined search direction)
    sigma_param = (1-alpha_p)^3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Rc = [Rcl, Rcq; Rce] and solve for the combined search direction from (28) and (27) with 
    % the above sigma (centering parameter)
    Rc = -x - sigma_param*mu_hat*g_dual(z, dimension_info);
    Rt = sigma_param*mu_hat/tau - kappa - dkappa_pred*dtau_pred/tau;
    h_hat = sparse(Rp_hat + A_hat*(H*Rd-Rc) + [sparse(m,1); Rt; 0]);
    % Find dy_hat = LHS_Schur_comp_eq \ h_hat; note that only h_hat has been updated
    % Use the same techniques as above
    Msp_inv_h_hat = Msp\h_hat;
    dy_hat = Msp_inv_h_hat - Msp_inv_U*(G\(U'*Msp_inv_h_hat));
    % Solve for dy_hat_actual for DEBUGGING
    % dy_hat = M_schur\h_hat;
    % disp(['error in dy_hat = ' num2str(norm(dy_hat-dy_hat_sm)/norm(dy_hat))]);
    % dy_hat = dy_hat_sm;
    
    dy = dy_hat(1:m);
    dtau = dy_hat(m+1);
    dtheta = dy_hat(m+2);
    dz = Rd - A_hat'*dy_hat; % 2nd equation od (27)
    dx = Rc - H*dz; % 3rd equation of (27)
    dkappa = Rt - (kappa/tau)*dtau - dkappa_pred*dtau_pred/tau; % 4th equation of (27)
    comb_dir = [dx; dy; dz; dtau; dkappa; dtheta];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Approximately find the step-length along the combined search direction and update the current iterate
    % alpha = (0.5+0.5*max(1-alpha_p,alpha_p))*find_alpha_max(x_bar, comb_dir, dimension_info);
    alpha = 0.8*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p, 0.8)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p,alpha_p)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = 0.98*find_alpha_max(x_bar, comb_dir, dimension_info);
    % x = x+alpha*dx; y = y + alpha*dy; z = z + alpha*dz; tau = tau + alpha*dtau; kappa = kappa + alpha*dkappa; theta = theta+alpha*dtheta;
    x_bar = x_bar + alpha*comb_dir;
    % Extract x, y, z, tau, theta, kappa
    x = x_bar(1:dim_x); y = x_bar(dim_x+1:dim_x+m); z = x_bar(dim_x+m+1:2*dim_x+m);  
    tau = x_bar(2*dim_x+m+1); kappa = x_bar(2*dim_x+m+2); theta = x_bar(2*dim_x+m+3);
    % Display some running information
    if it_count == 1
        disp(['size(A) = [' num2str(size(A,1)) ',' num2str(size(A,2)) '], density(A) = ' num2str(nnz(A)/(size(A,1)*size(A,2)))]);
        disp(['total_dim_l = ' num2str(Nl) ', total_dim_q = ' num2str(sum(Nq)) ', total_dim_e = ' num2str(3*Ne)]);
        %disp(['Schur_complement_equation_LHS size = [' num2str(size(LHS_Schur_comp_eq,1)) ', ' num2str(size(LHS_Schur_comp_eq,2)) ']']);
        %disp(['Initial density of it = ' num2str(nnz(LHS_Schur_comp_eq)/numel(LHS_Schur_comp_eq))]);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop started... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('  theta       sigma        dtheta        alpha         tau         kappa       iteration');
    end
    if rem(it_count,5) == 0
        fprintf('%10.4d | %10.4d | %10.4d | %10.4d | %10.4d | %10.4d | %4d \n', theta, sigma_param, full(dtheta), alpha, tau, kappa, it_count);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Termination conditions (from Section 5.4. in http://web.stanford.edu/~yyye/nonsymmhsdimp.pdf)
    % Check the 7 inequalities P, D, G, A, T, K, M
    bool_P = norm(A*x-tau*b,Inf) <= rel_eps*max(1,norm([A,b],Inf)); 
    bool_D = norm(A'*y+z-c*tau,Inf) <= rel_eps*max(1,norm([A',speye(dim_x),-c], Inf)); 
    bool_G = abs(-c'*x+b'*y-kappa) <= rel_eps*max(1, norm([-c',b',1],Inf)); 
    bool_A = abs(c'*x/tau - b'*y/tau) <= rel_eps*(1+abs(b'*y/tau)); 
    bool_T = tau <= rel_eps*(1e-2)*max(1,kappa); 
    bool_K = tau <= rel_eps*(1e-2)*min(1,kappa); 
    bool_M = theta <= rel_eps*(1e-2)*theta0; % bool_M = theta <= rel_eps*(1e-2)*theta0;
    if bool_P && bool_D && bool_A
        result_info.solution_status = 'optimal';
        display('The problem is primal and dual feasible.');
        display(['The approximate primal and dual optimal objectives are ' num2str(c'*x/tau) ' and ' num2str(b'*y/tau)]);
        break;
    elseif bool_P && bool_D && bool_G && bool_T
        is_primal_infeasible = (b'*y>0); is_dual_infeasible = (c'*x<0);
        if is_dual_infeasible && is_primal_infeasible
            result_info.solution_status = 'primal_and_dual_infeasible';
        elseif is_dual_infeasible
            result_info.solution_status = 'dual_infeasible';
        elseif is_primal_infeasible
            result_info.solution_status = 'primal_infeasible';
        end
        if is_dual_infeasible
            disp(['cTx = ' num2str(c'*x) ' < 0']);
            display('The problem is dual infeasible. See the info structure returned for a certificate.');
            % x_cert = x; c_vector = c; save('dual_infeas_cert.mat', 'x_cert', 'c_vector');
            idx_l = 1; idx_q = Nl+1; idx_e = Nl+sum(Nq)+1;
            for k = 1:size(blk,1)
                if blk{k,1} == 'l'
                    if indices_of_u(k) == 1
                        result_info.x_certificate_dual_infeasible{k} = x(idx_l:idx_l+blk{k,2}/2-1) - x(idx_l+blk{k,2}/2:idx_l+blk{k,2}-1);
                    else
                        result_info.x_certificate_dual_infeasible{k} = x(idx_l:idx_l+sum(blk{k,2})-1);
                    end
                    idx_l = idx_l + sum(blk{k,2});
                elseif blk{k,1} == 'q'
                    result_info.x_certificate_dual_infeasible{k} = x(idx_q:idx_q + sum(blk{k,2})-1);
                    idx_q = idx_q + sum(blk{k,2});
                elseif blk{k,1} == 'e'
                    result_info.x_certificate_dual_infeasible{k} = x(idx_e:idx_e + sum(blk{k,2})-1);
                    idx_e = idx_e + sum(blk{k,2});
                end
            end
        end
        if is_primal_infeasible
            disp(['bTy = ' num2str(b'*y) ' > 0']);
            display('The problem is primal infeasible. See the info structure returned for a certificate.');
            result_info.y_certificate_primal_infeasible = y;
        end
        break;
    elseif bool_K && bool_M
        display('The problem seems ill-posed and the solution returned might not make any sense!');
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check whether the program gets stuck halfway
    if dtheta*alpha == 0
        display('No more progress possible! The value of theta cannot decrease anymore! Algorithm terminates now!');
        res1 = norm(A*x/tau-b, Inf)/norm([A, b], Inf); res2 = norm(A'*y/tau + z/tau -c, Inf)/norm([A', speye(dim_x), c], Inf); res3 = abs(b'*y/tau - c'*x/tau)/(1+abs(b'*y/tau));
        display(['The relative residual norms (linear primal, linear dual, duality gap) are ' num2str(res1) ', ' num2str(res2) ' and ' num2str(res3)]);
        break;
    end
end

% After the main loop, check whether the maximum number of iterations has been reached
display(['Total number of iterations = ' num2str(it_count)]);
if it_count == max_iter_count
    display('Warning: maximum number of iterations reached without any conclusion!');
end

% Find the final solution (x_final, y_final, z_final)
x_final = x/tau; y_final = y/tau; z_final = z/tau;
% Construct the return tulpe
idx_l = 1; idx_q = Nl+1; idx_e = Nl+sum(Nq)+1;
for k = 1:size(blk,1)
    if blk{k,1} == 'l'
        if indices_of_u(k) == 1 
            % Need to construct x_return{k} differently if blk{k,1} was originally 'u'
            x_return{k} = x_final(idx_l:idx_l+blk{k,2}/2-1) - x_final(idx_l+blk{k,2}/2:idx_l+blk{k,2}-1);
            z_return{k} = zeros(blk{k,2}/2,1);
        else
            x_return{k} = x_final(idx_l:idx_l+blk{k,2}-1); 
            z_return{k} = z_final(idx_l:idx_l+blk{k,2}-1);
        end
        idx_l = idx_l + blk{k,2};
    elseif blk{k,1} == 'q'
        x_return{k} = x_final(idx_q:idx_q + sum(blk{k,2})-1);
        z_return{k} = z_final(idx_q:idx_q + sum(blk{k,2})-1); 
        idx_q = idx_q + sum(blk{k,2});
    elseif blk{k,1} == 'e'
        x_return{k} = x_final(idx_e:idx_e + sum(blk{k,2})-1);
        z_return{k} = z_final(idx_e:idx_e + sum(blk{k,2})-1); 
        idx_e = idx_e + sum(blk{k,2});
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% End of the main loop    
end
% Assign y_return
y_return = y_final;
obj_val = full([b'*y_final, c'*x_final]);
disp(['[dual_obj, primal_obj] = [' num2str(obj_val(1),5) ', ' num2str(obj_val(2),5) ']']);
% Construct the remaining entries of result_info
result_info.relative_duality_gap = abs(c'*x/tau - b'*y/tau)/(1+abs(b'*y/tau)); 
result_info.relative_primal_infeasibility = norm(A*x-tau*b,Inf)/max(1,norm([A,b],Inf));
result_info.relative_dual_infeasibility = norm(A'*y+z-c*tau,Inf)/max(1,norm([A',speye(dim_x),-c], Inf));
result_info.terminal_tau = tau;
result_info.terminal_kappa = kappa;
result_info.terminal_theta = theta;
% Print total running time
disp(['Total CPU time elapsed = ' num2str(cputime-t_begin) ' seconds']);
% End of function
format short;
end