function [obj_val, x_return,y_return,z_return, result_info] = hsd_lqeu_fast(blk, A_cell, c_cell, b, input_options)
% Calling syntax can be found here:
% https://github.com/gao-yuan-hangzhou/homogeneous_ipm_exp_cone/blob/master/README.md
format long;
warning('on', 'all');

% Include the path for subroutines
% All subroutines needed are in the folder "subroutines"
addpath('./subroutines');

disp('=== hsd_lqeu_fast (sparse LU with scaled partial pivoting) started... ===');
% This program makes use of [L,U,P,Q,R]=lu() 
% (sparse LU factorization with scaled partial pivoting)
% in solving the systems for the search directions.

% Get the current CPU time at the beginning
t_begin = cputime;

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
        Nq = [Nq; blk{k,2}]; % It's fine to let Nq change size in every iteration since size(blk,1) is usually small. 
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
    if blk{k,1} == 'l' cl = [cl; c_cell{k}];
    elseif blk{k,1} == 'q' cq = [cq; c_cell{k}];
    elseif blk{k,1} == 'e' ce = [ce; c_cell{k}];
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
% Possible initial values of tau, kappa and theta. Might not be the optimal setting.
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
% alpha_bar = (x'*z + tau*kappa)/theta;

% Formualte part of the coefficient matrix for the predictor search direction
% Note that the system for the predictor search direction can be compactly written as
% G_bar * dx_bar_p = R_bar_p, where G_bar = [G1; G2; G3] with
G1 = [-A, sparse(m,m), sparse(m,dim_x), b, sparse(m,1), -b_bar; 
    sparse(dim_x,dim_x), A', speye(dim_x), -c, sparse(dim_x,1), c_bar; c', -b', sparse(1,dim_x), 0, 1, -g_bar; 
    -c_bar', b_bar', sparse(1,dim_x), g_bar, 0, 0];
% and x_bar is the concatenated vector of all variables
x_bar = [x; y; z; tau; kappa; theta];

% Compute the dimension of the large x_bar for debugging
dim_x_bar = dim_x + m + dim_x + 3;

% Set the default relative accuracy and maximum iteration count
% if they are not set by the user
try 
    rel_eps = input_options.rel_eps;
catch err
    rel_eps = 1e-8;
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

% The main loop
for iter_idx =1:max_iter_count
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin of all termination conditions 
    % (from Section 5.4. in http://web.stanford.edu/~yyye/nonsymmhsdimp.pdf)
    % Extract x, y, z, tau, theta, kappa
    x = x_bar(1:dim_x); 
    y = x_bar(dim_x+1:dim_x+m); 
    z = x_bar(dim_x+m+1:2*dim_x+m);  
    tau = x_bar(2*dim_x+m+1); 
    kappa = x_bar(2*dim_x+m+2); 
    theta = x_bar(2*dim_x+m+3);
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
    % Form G_bar, R_bar_p and solve for the predictor search direction
    G2 = [sparse(1,2*dim_x+m), kappa, tau, 0]; 
    H = sparse(H_dual(z, dimension_info));
    % Assign H matrix to fix a Matlab memory allocation issue
    H = theta*H; 
    % Only the Hessian in G3 is updated, while the rest of G_bar remains unchanged throughout
    G3 = [speye(dim_x), sparse(dim_x, m), H, sparse(dim_x,3)]; 
    G_bar = [G1; G2; G3]; % DEBUGGING: save(['G_bar', num2str(cputime)], 'G_bar');
    R_bar_p = [sparse(m+dim_x+2,1); -tau*kappa; -x]; % R_bar_p = [sparse(m+dim_x+2,1); -tau*kappa; -x]; 
    % LU Factorization of G_bar using lu(G_bar, threshold)
    % Note that a smaller threshold gives less accurate LU factorization but 
    % is more memory-efficient
    try
        tbeginlu = cputime; 
        [L_lu, U_lu, P_lu, Q_lu, R_lu] = lu(G_bar, 0.01); 
        tlu = cputime - tbeginlu;
        % disp(['Time taken by the LU step = ' num2str(tlu)]);
    catch memory_error
        disp('Memory issue in factorizing G_bar...reduce the pivot threshold to 0.001');
        try 
            tbeginlu = cputime; 
            [L_lu, U_lu, P_lu, Q_lu, R_lu] = lu(G_bar, 0.001); 
            tlu = cputime - tbeginlu;
        catch memory_error
            disp('Memory issue in factorizing G_bar...reduce the pivot threshold to 0.0001');
            tbeginlu = cputime; 
            [L_lu, U_lu, P_lu, Q_lu, R_lu] = lu(G_bar, 0.0001); 
            tlu = cputime - tbeginlu;
        end
    end
    % DEBUGGING: LUmat = [L_lu, U_lu]; disp(['density of [L,U] = ' num2str(nnz(LUmat)/numel(LUmat))]);
    % DEBUGGING: disp(['nnz([L,U,P,Q,R])=' num2str(nnz(L_lu)+nnz(P_lu)+nnz(U_lu)+nnz(P_lu)+nnz(Q_lu)+nnz(R_lu))]);
    
    % Compute the predictor search direction by "solving" the system 
    % G_bar * pred_dir = R_bar_p, where
    % pred_dir = [dx_p, dy_p, dz_p, dtau_p, dkappa_p, dtheta_p]
    pred_dir = Q_lu*(U_lu\(L_lu\(P_lu*(R_lu\R_bar_p))));
    dtau_p = pred_dir(2*dim_x+m+1); 
    dkappa_p = pred_dir(2*dim_x+m+2);
    % DEBUGGING: dx_p = pred_dir(1:Nt); dy_p = pred_dir(Nt+1:Nt+m); dz_p = pred_dir(Nt+m+1: 2*Nt+m); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the step-length along the predictor search direction (approximately)
    alpha_p = find_alpha_max(x_bar, pred_dir, dimension_info);
    % Set sigma (weight of the linear combination of predidctor and corrector search directions)
    % A larger sigma gives a combined direction that is biased toward the 
    sigma = min(1,(1-alpha_p)^3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for the combined search direction
    R_bar = [sparse(dim_x+m+2, 1); -tau*kappa + sigma*theta - dtau_p*dkappa_p; -x - sigma*theta*g_dual(z, dimension_info)]; %R_bar = [sparse(dim_x+m+2, 1); -tau*kappa + sigma*theta - dtau_p*dkappa_p; -x - sigma*theta*g_dual(z, dimension_info)];
    comb_dir = Q_lu*(U_lu\(L_lu\(P_lu*(R_lu\R_bar)))); % comb_dir = G_bar\R_bar; % linsolve(G_bar,R_bar);
    % dx = comb_dir(1:Nt); dy = comb_dir(Nt+1:Nt+m); dz = comb_dir(Nt+m+1:2*Nt+m); dtau = comb_dir(2*Nt+m+1); dkappa = comb_dir(2*Nt+m+2); 
    dtheta = comb_dir(2*dim_x+m+3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Approximately find the step-length along the combined search direction and update the current iterate
    % alpha = (0.5+0.5*max(1-alpha_p,alpha_p))*find_alpha_max(x_bar, comb_dir, dimension_info);
    alpha_step = 0.8*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p, 0.8)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p,alpha_p)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = 0.98*find_alpha_max(x_bar, comb_dir, dimension_info);
    % x = x+alpha*dx; y = y + alpha*dy; z = z + alpha*dz; tau = tau + alpha*dtau; kappa = kappa + alpha*dkappa; theta = theta+alpha*dtheta;
    x_bar = x_bar + alpha_step*comb_dir;
    % disp([num2str(theta,5) ' | ' num2str(sigma,5) ' | ' num2str(dtheta,5) ' | ' num2str(alpha,5) ' | '  num2str(tau,5) ' | ' num2str(kappa,5) ' | ' num2str(nnz(G_bar)/numel(G_bar))]);
    if iter_idx == 1
        disp(['size(A) = [' num2str(size(A,1)) ',' num2str(size(A,2)) '], density(A) = ' num2str(nnz(A)/(size(A,1)*size(A,2)))]);
        disp(['total_dim_l = ' num2str(Nl) ', total_dim_q = ' num2str(sum(Nq)) ', total_dim_e = ' num2str(3*Ne)]);
        disp(['size(G_bar) = [' num2str(size(G_bar,1)) ', ' num2str(size(G_bar,2)) ']']);
        disp(['Initial nnz and density of G_bar = [' num2str(nnz(G_bar)) ', ' num2str(nnz(G_bar)/numel(G_bar)) ']']);
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop started... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('  theta       sigma        dtheta        alpha         tau         kappa       iteration');
    end
    if rem(iter_idx,5) == 0
        fprintf('%10.4d | %10.4d | %10.4d | %10.4d | %10.4d | %10.4d | %4d \n', theta, sigma, full(dtheta), alpha_step, tau, kappa, iter_idx);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check whether the program gets stuck halfway
    if false % dtheta*alpha_step == 0
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
disp(['Total CPU time elapsed = ' num2str(cputime-t_begin) ' seconds']);

% End of function
format short;
end
