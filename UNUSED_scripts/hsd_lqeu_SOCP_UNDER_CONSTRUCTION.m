function [obj_val, x_return,y_return,z_return, result_info] = hsd_lqeu(blk, A_cell, c_cell, b, input_options)
% Calling syntax can be found here:
% https://github.com/gao-yuan-hangzhou/homogeneous_ipm_exp_cone/blob/master/README.md
format long;

% Include the path for subroutines
% All subroutines needed are in the folder "subroutines"
warning('off', 'all');
addpath('./subroutines');
warning('on', 'all');

disp('=== hsd_lqeu_Schur_bicgstab (Schur complement equation solved using BiCGSTAB(l)) started... ===');
% This program makes use of [L,U,P,Q,R]=lu() 
% (sparse LU factorization with scaled partial pivoting)
% in solving the systems for the search directions.

% Get the current CPU time at the beginning
t_very_beginning = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether initial iterate is (partly) specified.
is_initial_x_given = false;
is_initial_y_given = false;
is_initial_z_given = false;

try 
    initial_x = input_options.initial_x;
    is_initial_x_given = true;
    disp('User specified initial_x!');
catch err
end

try 
    initial_y = input_options.initial_y;
    is_initial_y_given = true;
    disp('User specified initial_y!');
catch err
end

try initial_z = input_options.initial_z;
    is_initial_z_given = true;
    disp('User specified initial_z!');
catch err
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert all 'u' into 'l' using the simplest scheme x = x_pos - x_neg
% and keep track of the original indices of 'u'
% Also, if initial_x or initial_z is given, 
% convert their 'u' parts into 'l' parts
indices_of_u = zeros(size(blk,1),1);
for k = 1:size(blk,1)
    if blk{k,1} == 'u'
        indices_of_u(k) = 1;
        blk{k,1} = 'l'; blk{k,2} = 2*blk{k,2};
        A_cell{k} = [A_cell{k}, -A_cell{k}]; 
        c_cell{k} = [c_cell{k}; -c_cell{k}];
        if is_initial_x_given
            initial_x{k} = [max(initial_x{k},0); min(initial_x{k},0)]; 
        end
        if is_initial_z_given
            initial_z{k} = [max(initial_z{k},0); min(initial_z{k},0)];
        end
    end
end

% Say something about 'u'
if sum(indices_of_u)>0
    disp('Unrestricted blocks detected and are converted into linear blocks');
end

% Obtain the following dimensions:
% Nl (total dimension of linear part, after all 'u' converted into 'l'), 
% Nq = [q(1); ...; q(n_1)] (dimensions of second-order cones), 
% Ne (NUMBER of exponential cones)
% m (dimension of b, number of linear constraints, assuming rank(A)=m)
% Meanwhile, rearrange initial_x and initial_z into initial_x_vec and initial_z_vec
Nl = 0; Nq = []; Ne = 0; m = length(b);
initial_x_vec.l = []; initial_x_vec.q = []; initial_x_vec.e = [];
initial_z_vec.l = []; initial_z_vec.q = []; initial_z_vec.e = [];
for k = 1:size(blk,1)
    if blk{k,1} == 'l' 
        Nl = Nl+blk{k,2}; 
        if is_initial_x_given
            initial_x_vec.l = [initial_x_vec.l; initial_x{k}];
        end
        if is_initial_z_given
            initial_z_vec.l = [initial_z_vec.l; initial_z{k}];
        end
    elseif blk{k,1} == 'q' 
        Nq = [Nq; blk{k,2}];
        if is_initial_x_given
            initial_x_vec.q = [initial_x_vec.q; initial_x{k}];
        end
        if is_initial_z_given
            initial_z_vec.q = [initial_z_vec.q; initial_z{k}];
        end
    elseif blk{k,1} == 'e' 
        Ne = Ne + length(blk{k,2});
        if is_initial_x_given
            initial_x_vec.e = [initial_x_vec.e; initial_x{k}];
        end
        if is_initial_z_given
            initial_z_vec.e = [initial_z_vec.e; initial_z{k}];
        end
    else
        display('Error: blk{k,1} is not one of l, q, u or e!');
    end;
end;

initial_x_vec = [initial_x_vec.l; initial_x_vec.q; initial_x_vec.e];
initial_z_vec = [initial_z_vec.l; initial_z_vec.q; initial_z_vec.e];

% Compute the total dimension of x = [xl; xq; xe]
dim_x = Nl + sum(Nq) + 3*Ne;
n = dim_x; % alias
% Define structure dimension_info that contains Nl, Nq and Ne
dimension_info.l = Nl; dimension_info.q = Nq; dimension_info.e = Ne; dimension_info.m = m;
% Construct A = [Al, Aq, Ae]
Al = sparse(m,0); Aq = sparse(m,0); Ae = sparse(m,0);
for k=1:size(blk,1)
    if blk{k,1} == 'l' Al = [Al, A_cell{k}];
    elseif blk{k,1} == 'q' Aq = [Aq, A_cell{k}];
    elseif blk{k,1} == 'e' Ae = [Ae, A_cell{k}];
    else disp(['Error: blk{' num2str(k) ',1} is not one of l, q or e!']);
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

% Now initialize (x, y, z, tau, kappa, theta)
% If no initial values are given, 
% set x{j} = z{j} = conic identity (analytic centers) of K{j}
exp_cone_center = [-1.0151; 1.2590; 0.5560];
id_l = ones(Nl,1);
id_q = zeros(sum(Nq),1); 
for k = 1:length(Nq) 
    id_q(sum(Nq(1:k-1))+1) = 1; 
end;
id_e = repmat(exp_cone_center,Ne,1);
x0 = [id_l; id_q; id_e]; x = x0; y0 = zeros(m,1); y = y0; z0 = x0; z = z0;
% Juse some possible initial values of tau, kappa and theta
tau = 1; kappa = 1; theta0 = 1.0; theta = theta0;

% Set x = initial_x_vec, z = initial_z_vec, y = initial_y, if applicable
if is_initial_x_given
    x = initial_x_vec;
end
if is_initial_y_given
    y = initial_y;
end
if is_initial_z_given
    z = initial_z_vec;
end

% Compute the auxiliary parameters which completely define (HSD)
b_bar = (b*tau - A*x)/theta; 
c_bar = (c*tau - A'*y - z)/theta; 
g_bar = (c'*x - b'*y + kappa)/theta; 
alpha_bar = (x'*z + tau*kappa)/theta;

% Define the extented coefficient matrix
A_hat = [A; -c'; c_bar']; % A_hat is (m+2) by dim_x

% Extract dense columns from A_hat
try
    density_treshold = input_options.density_threshold;
catch err
    density_threshold = 0.08;  
end

% Find the dense columns in A_hat
density = @(X) nnz(X)/numel(X);
dense_col_bool = sparse(dim_x, 1);
A_hat_dc = sparse(m+2, dim_x);
A_hat_sp = A_hat;
for col_idx = 1:dim_x
    curr_col = A_hat(:, col_idx);
    if (density(curr_col) >= density_threshold)
        dense_col_bool(col_idx) = 1;
        A_hat_dc(:, col_idx) = A_hat(:, col_idx);
        A_hat_sp(:, col_idx) = 0;
    end
end
dense_col_indices = find(dense_col_bool);
num_dc = length(dense_col_indices);
disp(['Number of dense columns in A = ' num2str(num_dc)]);

% Find the permutation P such that A_hat * P = [A_hat_sp, A_hat_dc]
P_perm_mat = speye(dim_x, dim_x);
for k = 1:num_dc
    dc_idx = dense_col_indices(k);
    new_idx = n - num_dc + k;
    temp_col = P_perm_mat(:, dc_idx);
    P_perm_mat(:, dc_idx) = P_perm_mat(:, new_idx); 
    P_perm_mat(:, new_idx) = temp_col;
end

% Construct A_hat_sp and A_hat_dc
A_hat_sp = A_hat * P_perm_mat(:, 1:n-num_dc);
A_hat_dc = A_hat * P_perm_mat(:, n-num_dc+1:end);

% Set x_bar, the concatenated vector of all variables
x_bar = [x; y; z; tau; kappa; theta];

% Compute the dimension of the large x_bar (mainly for debugging)
dim_x_bar = 2*dim_x + m + 3;

% Set the default relative accuracy and maximum iteration count
% if they are not set by the user
try 
    rel_eps = input_options.rel_eps;
catch err
    rel_eps = 1e-6;
end

try
    max_iter_count = input_options.max_iter_count;
catch err
    max_iter_count = 500;
end

disp(['relative accuracy = ' num2str(rel_eps)]);
disp(['maximum iteration count = ' num2str(max_iter_count)]);

% Set the default solution status
result_info.solution_status = 'undetermined';

% Whether the systems for search directions are solved through LU factorization.
% The default is false (using Schur complement equation and )
is_lu = false;

% Set the BiCGSTAB tolerance and max_iter_bicgstab
tol_bicgstab = 1e-10;
max_iter_bicgstab = 8;

% Needed for decomposing the second-order cone Hessian
has_q = (sum(Nq)~=0);

% The main loop
for iter_idx =1:max_iter_count
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin of all termination conditions 
    % (from Section 5.4. in http://web.stanford.edu/~yyye/nonsymmhsdimp.pdf)
    % Extract x, y, z, tau, theta, kappa
    x = x_bar(1:dim_x); y = x_bar(dim_x+1:dim_x+m); z = x_bar(dim_x+m+1:2*dim_x+m);  
    tau = x_bar(2*dim_x+m+1); kappa = x_bar(2*dim_x+m+2); theta = x_bar(2*dim_x+m+3);
    % Check the 7 inequalities P, D, G, A, T, K, M
    % pinfeas = norm(A*x/tau-b,Inf)/max(1,norm([A,b],Inf));
    % dinfeas = norm(A'*y/tau+z/tau-c,Inf)/max(1,norm([A',speye(dim_x),-c], Inf));
    % reldualgap = abs(c'*x/tau - b'*y/tau)/(1+abs(b'*y/tau));
    bool_P = norm(A*x/tau-b,Inf) <= rel_eps*max(1,norm([A,b],Inf)); 
    bool_D = norm(A'*y/tau+z/tau-c,Inf) <= rel_eps*max(1, norm([A',speye(dim_x),-c], Inf)); 
    bool_G = abs(-c'*x+b'*y-kappa) <= rel_eps*max(1, norm([-c',b',1],Inf)); 
    bool_A = abs(c'*x/tau - b'*y/tau) <= rel_eps*(1+abs(b'*y/tau));
    bool_T = tau <= rel_eps*(1e-2)*max(1,kappa); 
    bool_K = tau <= rel_eps*(1e-2)*min(1,kappa); 
    bool_M = theta <= rel_eps*(1e-2)*theta0; % bool_M = theta <= rel_eps*(1e-2)*theta0;
    % disp (bool_P && bool_D && bool_A)
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
            display('The problem is dual infeasible. See info structure returned for a certificate.');
            % Construct and return a certificate of dual infeasibility
            idx_l = 1; 
            idx_q = Nl+1; 
            idx_e = Nl+sum(Nq)+1;
            for k = 1:size(blk,1)
                if blk{k,1} == 'l'
                    if indices_of_u(k) == 1
                        result_info.x_certificate_dual_infeasible{k} =...
                            x(idx_l:idx_l+blk{k,2}/2-1) - x(idx_l+blk{k,2}/2:idx_l+blk{k,2}-1);
                    else
                        result_info.x_certificate_dual_infeasible{k} =x(idx_l:idx_l+sum(blk{k,2})-1);
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
    % End of all termination conditions 
    % (except the case of numerical difficulty alpha_step * theta == 0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the search direction pred_dir = [dx; dy; dz; dtau; dkappa]
    % as in SDPT3 Guide: http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf
    % Construct the basic quantities A_hat, B_hat, the RHS components and so on
    y_hat = [y; tau; theta]; % y_hat has dimension (m+2)
    B_hat = [sparse(m,m), -b, b_bar; b', 0, g_bar; -b_bar', -g_bar, 0]; % B_hat is (m+2) by (m+2)
    mu_hat = theta; % mu_hat = (x'*z +tau*kappa)/(v_param+1) = theta; % mu_xz = x'*z/v;
    Rp_hat = [sparse(m,1); kappa; -alpha_bar] - A_hat*x - B_hat*y_hat;
    Rd = -A_hat'*y_hat - z;
    Rc = -x; % sigma = 0 for predictor direction system, as explained below
    Rt = - kappa; % Rt = mu_hat/tau - kappa;
    % The H matrix is the Hessian of the dual barrier evaluated at z
    % if iter_idx == 5 keyboard; end;
    % Compute the diagonal part of mu_hat*H_dual(z)
    % Definition: H = mu*H_dual(z) = E_mat + v*v'
    E_mat = E_dual(z, dimension_info); E_mat = mu_hat*E_mat; % to fix Matlab's memory allocation issue
    % Compute the rank-1 perturbation in the second-order cone Hessian
    if has_q
        v = sqrt(mu_hat)*v_dual(z, dimension_info); 
    else
        v = [];
    end
    % Define the linear operator H
    if has_q
        H = @(vec) E_mat*vec + v*(v'*vec);
    else
        H = @(vec) E_mat * vec;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now decompose M into M = Msp + U*D*U'; 
    % where Msp is sparse symmetric positive definite (nearly indefinite) and 
    % U*D*U' captures B_hat and the dense part of M due to dense columns of A_hat
    % First, formulate UB and DB such that UB*DB*UB' = B_hat
    U1 = [-b; 0; -g_bar/2]; U2 = [b_bar; g_bar/2; 0];
    V1 = [zeros(m,1); 1; 0]; V2 = [zeros(m+1,1); 1];
    UB = [U1, U2, V1, V2];
    DB = [zeros(2), eye(2); -eye(2), zeros(2)]; 
    DB_inv = [zeros(2), -eye(2); eye(2), zeros(2)];
    
    % Construct E_tilde_1 and E_tilde_2 such that 
    % blkdiag(E_tilde_1, E_tilde_2) = P'*E_mat*P
    E_tilde_1 = P_perm_mat(:, 1:n-num_dc)'*E_mat*P_perm_mat(:, 1:n-num_dc);
    E_tilde_2 = P_perm_mat(:, n-num_dc+1:end)'*E_mat*P_perm_mat(:, n-num_dc+1:end);

    % Set M, Msp, U, D such that 
    % M = Msp + U*D*U' 
    Msp = A_hat_sp * E_tilde_1 * A_hat_sp' + diag([sparse(m,1); kappa/tau; 0]);
    Msp_pert = Msp+1e-15*norm(Msp, Inf)*speye(m+2); % chol(Msp_pert)
    
    if has_q
        U = [UB, A_hat_dc, A_hat*v];
        D = blkdiag(DB, E_tilde_2, 1);
        D_inv = blkdiag(DB_inv, inv(E_tilde_2), 1);
    else
        U = [UB, A_hat_dc];
        D = blkdiag(DB, E_tilde_2);
        D_inv = blkdiag(DB_inv, inv(E_tilde_2));
    end
    
    % Define M = Msp + U*D*U' as a linear map since computing U*D*U' is time-consuming
    M = @(yy) Msp*yy + U*D*(U'*yy);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve equation (28) and then (27) with sigma = 0 for the predictor direction
    % (28) is M*y_hat = h_hat (with sigma = 0)
    % Construct the right hand side
    h_hat = sparse(Rp_hat + A_hat*(H(Rd)-Rc) + [sparse(m,1); Rt; 0]); 
    
    % Solve for dy_hat using Sherman-Morrison formula
    % Calculate Msp\h_hat = R\(R'\h_hat), Msp\U
    % and the 4-by-4 matrix G = D_inv + U'*inv(Msp)*U
    if ~is_lu
        try
            Lchol = chol(Msp_pert);
            Msp_inv_chol = @(X) Lchol\(Lchol'\X);
            Msp_inv_U = Msp_inv_chol(U);
            % Msp_inv_h_hat = Msp_inv_chol(h_hat);
            Gmat = D_inv + U'*Msp_inv_U;
            [LG, UG] = lu(Gmat);
            % Implicity calculate P = eye(m+2) - inv(Msp)*U*G_inv*U' and then y_hat = P*w_hat
            % P = speye(m+2) - Msp_inv_U*(G\U');
            % Calculate dy_hat_pred = LHS_Schur_comp_eq \ h_hat; % (28)
            % Mathematically, dy_hat_pred = Msp_inv_h_hat - Msp_inv_U*(G\(U'*Msp_inv_h_hat));
            precond_func = @(rhs) precond_M_inv(rhs, Lchol, U, LG, UG, Msp_inv_U);
            [dy_hat_pred, flag] = bicgstab(M, h_hat, tol_bicgstab, max_iter_bicgstab, precond_func);
            if flag ~= 0
                throw(exception);
            end
        catch err
            disp(['Preconditioned BiCGSTAB fails! Switch to LU factorization of M_aug at iteration ', num2str(iter_idx)]);
            keyboard;
            is_lu = true;
        end
    end
    
    % If is_lu = true, solve equation (21) in SDPT3 guide  
    if is_lu
        h_hat_aug = [h_hat; zeros(num_dc+4+has_q,1)]; % D is (s+4)-by-(s+4), where s is the number of dense columns
        M_aug = [Msp_pert, U; U', -D_inv]; % M_aug is (m+2+s+4)-by-(m+2+s+4), slightly larger than Msp
        [L_lu, U_lu, P_lu, Q_lu, R_lu] = lu(M_aug, 0.05);
        dy_hat_aug = Q_lu*(U_lu\(L_lu\(P_lu*(R_lu\h_hat_aug))));
        dy_hat_pred = dy_hat_aug(1:m+2);
    end
    
    % Compute dy_pred, dx_pred, ...
    dy_pred = dy_hat_pred(1:m);
    dtau_pred = dy_hat_pred(m+1);
    dtheta_pred = dy_hat_pred(m+2);
    dz_pred = Rd - A_hat'*dy_hat_pred; % 2nd equation of (27)
    dx_pred = Rc - H(dz_pred); % 3rd equation of (27)
    dkappa_pred = Rt - (kappa/tau)*dtau_pred; % 4th equation of (27)
    pred_dir = [dx_pred; dy_pred; dz_pred; dtau_pred; dkappa_pred; dtheta_pred];

    % if (iter_idx==5) keyboard; end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Given current iterate x_bar and pred_dir, approximately compute alpha_p
    x_bar = [x; y; z; tau; kappa; theta];
    alpha_p = find_alpha_max(x_bar, pred_dir, dimension_info);
    % Set sigma (parameter for the system for the combined search direction)
    sigma = min(1,(1-alpha_p)^3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find Rc = [Rcl, Rcq; Rce] and solve for the combined search direction from (28) and (27) with 
    % the above sigma (centering parameter)
    % Solve for the combined direction using either BiCGSTAB or LU factorization
    Rc = -x - sigma*mu_hat*g_dual(z, dimension_info);
    Rt = sigma*mu_hat/tau - kappa - dkappa_pred*dtau_pred/tau;
    h_hat = sparse(Rp_hat + A_hat*(H(Rd)-Rc) + [sparse(m,1); Rt; 0]);
    if is_lu
        h_hat_aug = [h_hat; zeros(num_dc+4+has_q,1)];
        dy_hat_aug = Q_lu*(U_lu\(L_lu\(P_lu*(R_lu\h_hat_aug))));
        dy_hat = dy_hat_aug(1:m+2);
    else % Otherwise, still use preconditioned BiGSTAB, Cholesky and LU factors already computed
        % dy_hat = bicgstab(M, h_hat, tol_bicgstab, max_iter_bicgstab, precond_func);
        [dy_hat, flag] = bicgstab(M, h_hat, tol_bicgstab, max_iter_bicgstab, precond_func);
    end      
    dy = dy_hat(1:m);
    dtau = dy_hat(m+1);
    dtheta = dy_hat(m+2);
    dz = Rd - A_hat'*dy_hat; % 2nd equation od (27)
    dx = Rc - H(dz); % 3rd equation of (27)
    dkappa = Rt - (kappa/tau)*dtau - dkappa_pred*dtau_pred/tau; % 4th equation of (27)
    comb_dir = [dx; dy; dz; dtau; dkappa; dtheta];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Approximately find the step-length along the combined search direction and update the current iterate
    % alpha = (0.5+0.5*max(1-alpha_p,alpha_p))*find_alpha_max(x_bar, comb_dir, dimension_info);
    alpha_combined = 0.8*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p, 0.8)*find_alpha_max(x_bar, comb_dir, dimension_info);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the iterate
    % x = x+alpha*dx; y = y + alpha*dy; z = z + alpha*dz; tau = tau + alpha*dtau; kappa = kappa + alpha*dkappa; theta = theta+alpha*dtheta;
    x_bar = x_bar + alpha_combined*comb_dir;
    % disp([num2str(theta,5) ' | ' num2str(sigma,5) ' | ' num2str(dtheta,5) ' | ' num2str(alpha,5) ' | '  num2str(tau,5) ' | ' num2str(kappa,5) ' | ' num2str(nnz(G_bar)/numel(G_bar))]);
    if iter_idx == 1
        disp(['size(A) = [' num2str(size(A,1)) ',' num2str(size(A,2)) '], density(A) = ' num2str(nnz(A)/(size(A,1)*size(A,2)))]);
        disp(['total_dim_l = ' num2str(Nl) ', total_dim_q = ' num2str(sum(Nq)) ', total_dim_e = ' num2str(3*Ne)]);
        % disp(['size(G_bar) = [' num2str(dim_x_bar) ', ' num2str(dim_x_bar) ']']);
        %disp(['Initial nnz and density of G_bar = [' num2str(nnz(G_bar)) ', ' num2str(nnz(G_bar)/numel(G_bar)) ']']);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp(' theta  pinfeas dinfeas dualgap  primalobj    dualobj     steplen   tau    kappa  iter');
    end
    if rem(iter_idx,5) == 0
        pinfeas = norm(A*x/tau-b,Inf)/max(1,norm([A,b],Inf)); % primal infeasibility
        dinfeas = norm(A'*y/tau+z/tau-c,Inf)/max(1,norm([A',speye(dim_x),-c], Inf)); % dual infeasibility
        reldualgap = abs(c'*x/tau - b'*y/tau)/(1+abs(b'*y/tau)); % relative duality gap
        fprintf('%5.1d|%5.1d|%5.1d|%5.1d|%12.6d|%12.6d|%5.1d|%5.1d|%5.1d|%3d \n', ...
            theta, pinfeas, dinfeas, reldualgap, c'*x/tau, b'*y/tau, alpha_combined, tau, kappa, iter_idx);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check whether the program gets stuck halfway
    if dtheta*alpha_combined == 0
        display('No more progress possible! The value of theta cannot decrease anymore! Algorithm terminates now!');
        res1 = norm(A*x/tau-b, Inf)/norm([A, b], Inf); 
        res2 = norm(A'*y/tau + z/tau -c, Inf)/norm([A' speye(dim_x) c], Inf); 
        res3 = abs(b'*y/tau - c'*x/tau)/(1+abs(b'*y/tau));
        display(['The relative residual norms (linear primal, linear dual, duality gap) are ' ...
            num2str(res1) ', ' num2str(res2) ' and ' num2str(res3)]);
        break;
    end
end

% After the main loop, check whether the maximum number of iterations has been reached
% Total number of iterations = terminal iteration count - 1
% since termination conditions go before the main body of an iteration
actual_iter_count = iter_idx - 1;
display(['Total number of iterations = ' num2str(actual_iter_count)]);
if iter_idx == max_iter_count
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
% End of the main loop    
end

% Assign y_return
y_return = y_final;
% Take the dual objective value as the approximate optimal objective value
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
disp(['Total CPU time elapsed = ' num2str(cputime-t_very_beginning) ' seconds']);

% End of function
format short;
end
