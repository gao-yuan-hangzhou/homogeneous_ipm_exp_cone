function [obj_val, x_return,y_return,z_return, result_info] = hsd_lqeu(blk, At, c_input, b, rel_eps, max_iter_count)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves problems of the following form
% min c'x s.t. Ax + (Au)(xu)=b, x in K, where
% K = (R_+)^(n_l) × Q(q(1))×Q(q(2))×...Q(q(n_q)) × (K_exp)^(n_e)
% R_+ = {x: x>=0}.
% Q(n) = {x in R^n: x(1)>=||x(2:n)||}, where || || is the Euclidean 2-norm.
% K_exp = closure{(x1,x2,x3): x2>=0, x3>0, exp(x1/x3)<=x2/x3} is the exponential cone.
% The variable xu is free (unconstrained).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that we support input of the following form:
% min sum_j (c(j)'x(j)) s.t. sum_j (A(j)x(j)) = b, x(j) in K(j) (or free), j=1,2,...,N
% where K(j) can be a nonnegative orthant of arbitrary dimension, a Lorentz cone, product of Lorentz
% cones, product of exponential cones or just R^n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The exact input format is described as follows.
% [obj,X,y,Z,info] = hsd_lqe(blk,At,c,b,rel_eps)
% Input arguments:
% blk: a cell array describing the dimensions of input data
% At: the data for A in Ax=b
% b and c: as defined above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: Suppose one has 
% x(1) in (R_+)^15, x(2) in Q(3)×Q(4)×Q(6), x(3) in (K_exp)^2, x(4) in Q(8), 
% x(5) in (K_exp)^7, x(6) in R^9, m=23.
% Then necessarily c = (c{1}, c{2}, c{3}, c{4}, c{5}) is a vector of dimension 
% 15 + (3+4+6) + 3*2 + 8 + 3*7 + 9 = 72
% b is a vector of dimension m=23.
% The input arguments are
% blk{1,1} = 'l'; blk{1,2} = 15;        At{1} is a sparse 23-by-15 matrix
% blk{2,1} = 'q'; blk{2,2} = [3; 4; 6]; At{2} is a sparse 23-by-13 matrix  (3+4+6=13)
% blk{3,1} = 'e'; blk{3,2} = 2;         At{3} is a sparse matrix of dimension 23-by-6
% blk{4,1} = 'q'; blk{4,2} = [8];       At{4} is a sparse matrix of dimension 23-by-8
% blk{5,1} = 'e'; blk{5,2} = 7;         At{5} is a sparse matrix of dimension 23-by-21
% blk{6,1} = 'u'; blk{6,2} = 9;         At{6} is a sparse 23-by-9 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First of all, convert all 'u' into 'l' using the simplest scheme x = x_pos - x_neg
indices_of_u = zeros(size(blk,1),1);
for k = 1:size(blk,1)
    if blk{k,1} == 'u'
        indices_of_u(k) = 1;
        blk{k,1} = 'l'; blk{k,2} = 2*blk{k,2};
        At{k} = [At{k}, -At{k}]; 
        c_input{k} = [c_input{k}; -c_input{k}];
    end
end

% Obtain Nl, Nq = [q(1), ..., q(n_1)], Ne and m
Nl = 0; Nq = []; Ne = 0; m = length(b);
for k = 1:size(blk,1)
    if blk{k,1} == 'l' Nl = Nl+blk{k,2}; 
    elseif blk{k,1} == 'q' Nq = [Nq; blk{k,2}];
    elseif blk{k,1} == 'e' Ne = Ne + blk{k,2};
    else display('Error: blk{k,1} is not one of l, q or e!');
    end;
end;
% Define Nt, the total dimension of x = [xl; xq; xe]
Nt = Nl + sum(Nq) + 3*Ne;
% Define structure dimension_info that contains Nl, Nq and Ne
dimension_info.l = Nl; dimension_info.q = Nq; dimension_info.e = Ne; dimension_info.m = m;
% Construct A = [Al, Aq, Ae]
Al = sparse(m,0); Aq = sparse(m,0); Ae = sparse(m,0);
for k=1:size(blk,1)
    if blk{k,1} == 'l' Al = [Al, At{k}];
    elseif blk{k,1} == 'q' Aq = [Aq, At{k}];
    elseif blk{k,1} == 'e' Ae = [Ae, At{k}];
    else disp(['Error: blk{' num2str(k) ',1} is not one of l, q or e!']);
    end;
end;
A = [Al, Aq, Ae];

% Contruct c = [cl; cq; ce]
cl = sparse(0,1); cq = sparse(0,1); ce = sparse(0,1);
for k=1:size(blk,1)
    if blk{k,1} == 'l' cl = [cl; c_input{k}];
    elseif blk{k,1} == 'q' cq = [cq; c_input{k}];
    elseif blk{k,1} == 'e' ce = [ce; c_input{k}];
    else display('Error: blk{k,1} is not one of l, q or e!');
    end;
end;
c = [cl; cq; ce];

% Check whether Ax=b has a solution
[L_A, U_A] = lu(A); [L_full,U_full] = lu([A,b]);
if sprank(U_A) < sprank(U_full)
    display('Error: Ax=b has no solution! The primal problem is infeasible.'); return;
elseif sprank(U_A) < m
    display('Error: A is not FULL ROW RANK!');
    return;
end

% Now initialize (x, y, z, tau, kappa, theta)
iota = [-1.0151; 1.2590; 0.5560];
id_l = ones(Nl,1);
id_q = zeros(sum(Nq),1); for k = 1:length(Nq) id_q(sum(Nq(1:k-1))+1) = 1; end;
id_e = repmat(iota,Ne,1);
x0 = [id_l; id_q; id_e]; x = x0; y0 = zeros(m,1); y = y0; z0 = x0; z = z0;
tau = 1; kappa = 1; theta0 = 1; theta = theta0;
% Compute the auxiliary parameters which completely define (HSD)
b_bar = (b*tau - A*x)/theta; 
c_bar = (c*tau - A'*y - z)/theta; 
g_bar = (c'*x - b'*y + kappa)/theta; 
alpha_bar = (x'*z + tau*kappa)/theta;

% Formualte part of the coefficient matrix for the predictor search direction
% Note that the system for the predictor search direction can be compactly written as
% G_bar * dx_bar_p = R_bar_p, where G_bar = [G1; G2; G3] with
G1 = [-A, sparse(m,m), sparse(m,Nt), b, sparse(m,1), -b_bar; sparse(Nt,Nt), A', speye(Nt), -c, sparse(Nt,1), c_bar; c', -b', zeros(1,Nt), 0, 1, -g_bar; -c_bar', b_bar', sparse(1,Nt), g_bar, 0, 0];
% and x_bar is the concatenated vector of all variables
x_bar = [x; y; z; tau; kappa; theta];

% Set the default relative accuracy
if nargin < 5
    rel_eps = 1e-8;
end

result_info.solution_status = 'undetermined';

if nargin < 6
    max_iter_count = 500;
end

disp(['Algorithm begins with size(A) = [' num2str(size(A,1)) ',' num2str(size(A,2)) '], density(A) = ' num2str(nnz(A)/(size(A,1)*size(A,2)))]);
disp(['Nl = ' num2str(Nl) ', sum(Nq) = ' num2str(sum(Nq)) ', 3*Ne = ' num2str(3*Ne)]);
disp('theta | sigma | dtheta | alpha | tau | kappa | density(G_bar)');
for it_cout =1:max_iter_count
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Form G_bar, R_bar_p and solve for the predictor search direction
    G2 = [sparse(1,2*Nt+m), kappa, tau, 0]; G3 = [speye(Nt), sparse(Nt, m), theta*H_dual(z, dimension_info), sparse(Nt,3)];
    G_bar = [G1; G2; G3];
    R_bar_p = [sparse(m+Nt+2,1); -tau*kappa; -x];
    pred_dir = G_bar\R_bar_p; %linsolve(G_bar,R_bar_p);
    % dx_p = pred_dir(1:Nt); dy_p = pred_dir(Nt+1:Nt+m); dz_p = pred_dir(Nt+m+1: 2*Nt+m); 
    dtau_p = pred_dir(2*Nt+m+1); dkappa_p = pred_dir(2*Nt+m+2); % dtheta_p = pred_dir(2*Nt+m+3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute alpha_p (approximately)
    alpha_p = find_alpha_max(x_bar, pred_dir, dimension_info);
    % Set sigma (parameter for the system for the combined search direction)
    sigma = min(1,(1-alpha_p)^3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for the combined search direction
    R_bar = [sparse(Nt+m+2, 1); -tau*kappa + sigma*theta - dtau_p*dkappa_p; -x - sigma*theta*g_dual(z, dimension_info)];
    comb_dir = G_bar\R_bar; % linsolve(G_bar,R_bar);
    % dx = comb_dir(1:Nt); dy = comb_dir(Nt+1:Nt+m); dz = comb_dir(Nt+m+1:2*Nt+m); dtau = comb_dir(2*Nt+m+1); dkappa = comb_dir(2*Nt+m+2); 
    dtheta = comb_dir(2*Nt+m+3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Approximately find the step-length along the combined search direction and update the current iterate
    % alpha = (0.5+0.5*max(1-alpha_p,alpha_p))*find_alpha_max(x_bar, comb_dir, dimension_info);
    alpha = 0.8*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p, 0.8)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = max(1-alpha_p,alpha_p)*find_alpha_max(x_bar, comb_dir, dimension_info);
    % alpha = 0.98*find_alpha_max(x_bar, comb_dir, dimension_info);
    % x = x+alpha*dx; y = y + alpha*dy; z = z + alpha*dz; tau = tau + alpha*dtau; kappa = kappa + alpha*dkappa; theta = theta+alpha*dtheta;
    x_bar = x_bar + alpha*comb_dir;
    display([num2str(theta,5) ' | ' num2str(sigma,5) ' | ' num2str(dtheta,5) ' | ' num2str(alpha,5) ' | '  num2str(tau,5) ' | ' num2str(kappa,5) ' | ' num2str(nnz(G_bar)/numel(G_bar))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Termination conditions (from Section 5.4. in http://web.stanford.edu/~yyye/nonsymmhsdimp.pdf)
    % Extract x, y, z, tau, theta, kappa
    x = x_bar(1:Nt); y = x_bar(Nt+1:Nt+m); z = x_bar(Nt+m+1:2*Nt+m);  
    tau = x_bar(2*Nt+m+1); kappa = x_bar(2*Nt+m+2); theta = x_bar(2*Nt+m+3);
    % Check the 7 inequalities P, D, G, A, T, K, M
    bool_P = norm(A*x-tau*b,Inf) <= rel_eps*max(1,norm([A,b],Inf)); 
    bool_D = norm(A'*y+z-c*tau,Inf) <= rel_eps*max(1,norm([A',speye(Nt),-c], Inf)); 
    bool_G = abs(-c'*x+b'*y-kappa) <= rel_eps*max(1, norm([-c',b',1],Inf)); 
    bool_A = abs(c'*x/tau - b'*y/tau) <= rel_eps*(1+abs(b'*y/tau)); 
    bool_T = tau <= rel_eps*(1e-2)*max(1,kappa); 
    bool_K = tau <= rel_eps*(1e-2)*min(1,kappa); 
    bool_M = theta <= rel_eps*(1e-2)*theta0;
    if bool_P && bool_D && bool_A
        result_info.solution_status = 'optimal';
        display('The problem is primal and dual feasible.');
        display(['The approximate primal and dual optimal objectives are ' num2str(c'*x/tau) ' and ' num2str(b'*y/tau)]);
        break;
    elseif bool_P && bool_D && bool_G && bool_T
        if c'*x < 0
            result_info.solution_status = 'dual infeasible';
            display('The problem is dual infeasible. See the info structure returned for a certificate.');
            idx_l = 1; idx_q = Nl+1; idx_e = Nl+sum(Nq)+1;
            for k = 1:size(blk,1)
                if blk{k,1} == 'l'
                    result_info.x_certificate_dual_infeasible{k} = x(idx_l:idx_l+blk{k,2}-1);
                    idx_l = idx_l + blk{k,2};
                elseif blk{k,1} == 'q'
                    result_info.x_certificate_dual_infeasible{k} = x(idx_q:idx_q + sum(blk{k,2})-1);
                    idx_q = idx_q + sum(blk{k,2});
                elseif blk{k,1} == 'e'
                    result_info.x_certificate_dual_infeasible{k} = x(idx_e:idx_e + 3*blk{k,2}-1);
                    idx_e = idx_e + 3*blk{k,2};
                end
            end
        elseif b'*y > 0
            result_info.solution_status = 'primal infeasible';
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
    if abs(dtheta*alpha)/theta <= 1e-10
        display('No more progress possible! The value of theta cannot decrease anymore! Algorithm terminates now!');
        res1 = norm(A*x/tau-b, Inf)/norm([A, b], Inf); res2 = norm(A'*y/tau + z/tau -c, Inf)/norm([A' eye(Nt) c], Inf); res3 = abs(b'*y/tau - c'*x/tau)/(1+abs(b'*y/tau));
        display(['The relative residual norms (linear primal, linear dual, duality gap) are ' num2str(res1) ', ' num2str(res2) ' and ' num2str(res3)]);
        break;
    end
end

% After the main loop, check whether the maximum number of iterations has been reached
display(['Total number of iterations = ' num2str(it_cout)]);
if it_cout == max_iter_count
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
        x_return{k} = x_final(idx_e:idx_e + 3*blk{k,2}-1);
        z_return{k} = z_final(idx_e:idx_e + 3*blk{k,2}-1); 
        idx_e = idx_e + 3*blk{k,2};
    end
% End of the main loop    
end
% Assign y_return
y_return = y_final;
% Take the dual objective value as the approximate optimal objective value
obj_val = [c'*x_final, b'*y_final];
% Construct the remaining entries of result_info
result_info.relative_duality_gap = abs(c'*x/tau - b'*y/tau)/(1+abs(b'*y/tau)); 
result_info.relative_primal_infeasibility = norm(A*x-tau*b,Inf)/max(1,norm([A,b],Inf));
result_info.relative_dual_infeasibility = norm(A'*y+z-c*tau,Inf)/max(1,norm([A',eye(Nt),-c], Inf));
result_info.terminal_tau = tau;
result_info.terminal_kappa = kappa;
result_info.terminal_theta = theta;
% End of function
end
