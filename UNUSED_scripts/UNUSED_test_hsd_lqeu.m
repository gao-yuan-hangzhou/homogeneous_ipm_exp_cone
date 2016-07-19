% We generate (feasible) instances of conic programs with 
% linear, second-order and exponential cone constraints
% to test hsd_lqe.m

% Choose whether to assure feasibility of the generated instance
is_cheating = true;
% Example: Suppose one has 
% x(1) in (R_+)^15, 
% x(2) in Q(3)×Q(4)×Q(6), 
% x(3) in (K_exp)^2, 
% x(4) in Q(8), 
% x(5) in (K_exp)^7,
% x(6) in (R_+)^9,
% m=23.
% Then necessarily,
% c = (c{1}, c{2}, c{3}, c{4}, c{5})
% is a vector of dimension 15 + (3+4+6) + 3*2 + 8 + 3*7 + 9 = 72
% b is a vector of dimension m=23.
% The input arguments are
% blk{1,1} = 'l'; blk{1,2} = 15; At{1} is a sparse 23-by-15 matrix
% blk{2,1} = 'q'; blk{2,2} = [3; 4; 6]; At{2} is a sparse 23-by-13 matrix  (3+4+6=13)
% blk{3,1} = 'e'; blk{3,2} = 2; At{3} is a sparse matrix of dimension 23-by-6
% blk{4,1} = 'q'; blk{4,2} = [8]; At{4} is a sparse matrix of dimension 23-by-8
% blk{5,1} = 'e'; blk{5,2} = 7; At{5} is a sparse matrix of dimension 23-by-21
% blk{6,1} = 'l'; blk{6,2} = 9; At{1} is a sparse 23-by-9 matrix

% Construct random blk, At, c and b
clear blk; clear c_input; clear At; m = 137;
blk{1,1} = 'l'; blk{1,2} = 1; At{1} = sprandn(m,sum(blk{1,2}),0.05); c_input{1} = randn(blk{1,2},1);
blk{2,1} = 'q'; blk{2,2} = [3;4;6;7;2;2;5;2;5;4;3;7;6;5;4;5;7;8;9;5;10]; At{2} = sprandn(m, sum(blk{2,2}), 0.1); c_input{2} = randn(sum(blk{2,2}),1);
blk{3,1} = 'e'; blk{3,2} = 100; At{3} = sprandn(m, 3*blk{3,2}, 0.15); c_input{3} = randn(3*blk{3,2}, 1);
blk{4,1} = 'q'; blk{4,2} = [8; 11; 6]; At{4} = sprandn(m, sum(blk{4,2}), 0.12); c_input{4} = randn(sum(blk{4,2}), 1);
blk{5,1} = 'l'; blk{5,2} = 1; At{5} = sprandn(m, blk{5,2}, 0.08); c_input{5} = randn(blk{5,2}, 1);
blk{6,1} = 'l'; blk{6,2} = 1; At{6} = sprandn(m, blk{6,2}, 0.07); c_input{6} = randn(blk{6,2}, 1);
blk{7,1} = 'e'; blk{7,2} = 140; At{7} = sprandn(m, 3*blk{7,2}, 0.11); c_input{7} = randn(3*blk{7,2},1);
      
b = randn(m,1);

% Choose whether to have a feasible instance
if is_cheating
    % Obtain the dimensions of different cones
    Nl = 0; Nq = []; Ne = 0; m = length(b);
    for k = 1:size(blk,1) 
        if blk{k,1} == 'l' Nl = Nl+blk{k,2}; 
        elseif blk{k,1} == 'q' Nq = [Nq; blk{k,2}];
        elseif blk{k,1} == 'e' Ne = Ne + blk{k,2};
        elseif blk{k,1} == 'u' Nl = Nl + 2*blk{k,2};
        end;
    end;
    Ntotal = Nl + sum(Nq) + 3*Ne;
    display(['Constructing a feasible instance...']);
    display(['Total dimension of x(linear) = ' num2str(Nl)]);
    display(['Total dimension of x(quadratic) = ' num2str(sum(Nq))]);
    display(['Total dimension of x(exponential) = ' num2str(3*Ne)]);
    display(['Number of equality constraints = ' num2str(m)]);
    % Construct [A, c, b] such that Ax=b, A'y+z=c is primal and dual feasible
    x_exp_primal_interior = [1;5;3]; z_exp_dual_interior = [-5;2;1]; b = zeros(m,1); y = rand(m,1);
    for k = 1:size(blk,1)
        if blk{k,1} == 'l'  
            b = b + At{k}*exp(rand(blk{k,2},1)); c_input{k} = At{k}'*y + exp(rand(blk{k,2},1));
        elseif blk{k,1} == 'q' 
            b = b + At{k}*exp(rand())*get_id_soc(blk{k,2}); c_input{k} = At{k}'*y + exp(rand())*get_id_soc(blk{k,2});
        elseif blk{k,1} == 'e' 
            b = b + At{k}*repmat(x_exp_primal_interior, [blk{k,2},1]); c_input{k} = At{k}'*y + repmat(z_exp_dual_interior, [blk{k,2},1]);
        end
    end
end

[obj_val, x_return,y_return,z_return, info] = hsd_lqeu(blk, At, c_input, b, 1e-8, 500);
display(obj_val);

% %%%%%%%%%% Function goes here %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Obtain Nl, Nq = [q(1), ..., q(n_1)], Ne and m
% Nl = 0; Nq = []; Ne = 0; m = length(b);
% for k = 1:size(blk,1) 
%     if blk{k,1} == 'l' Nl = Nl+blk{k,2}; 
%     elseif blk{k,1} == 'q' Nq = [Nq; blk{k,2}];
%     elseif blk{k,1} == 'e' Ne = Ne + blk{k,2};
%     else display('Error: blk{k,1} is not one of l, q or e!');
%     end;
% end;
% % Define Nt, the total dimension of x = [xl; xq; xe]
% Nt = Nl + sum(Nq) + 3*Ne;
% % Define structure dimension_info that contains Nl, Nq and Ne
% dimension_info.l = Nl; dimension_info.q = Nq; dimension_info.e = Ne; dimension_info.m = m;
% % Construct A = [Al, Aq, Ae]
% Al = sparse(m,0); Aq = sparse(m,0); Ae = sparse(m,0);
% for k=1:size(blk,1)
%     if blk{k,1} == 'l' Al = [Al, At{k}];
%     elseif blk{k,1} == 'q' Aq = [Aq, At{k}];
%     elseif blk{k,1} == 'e' Ae = [Ae, At{k}];
%     else display('Error: blk{k,1} is not one of l, q or e!');
%     end;
% end;
% A = [Al, Aq, Ae];
% 
% % Contruct c = [cl; cq; ce]
% cl = sparse(0,1); cq = sparse(0,1); ce = sparse(0,1);
% for k=1:size(blk,1)
%     if blk{k,1} == 'l' cl = [cl; c{k}];
%     elseif blk{k,1} == 'q' cq = [cq; c{k}];
%     elseif blk{k,1} == 'e' ce = [ce; c{k}];
%     else display('Error: blk{k,1} is not one of l, q or e!');
%     end;
% end;
% c = [cl; cq; ce];
% 
% % Check whether Ax=b has a solution
% [L_A, U_A] = lu(A); [L_full,U_full] = lu([A,b]);
% if sprank(U_A) < sprank(U_full)
%     display('Warning: Ax=b has no solution! The problem is primal infeasible.'); return;
% elseif sprank(U_A) < m
%     display('Warning: A is not FULL ROW RANK!');
% end
% 
% % Now initialize (x, y, z, tau, kappa, theta)
% iota = [-1.0151; 1.2590; 0.5560];
% id_l = ones(Nl,1);
% id_q = zeros(sum(Nq),1); for k = 1:length(Nq) id_q(sum(Nq(1:k-1))+1) = 1; end;
% id_e = repmat(iota,Ne,1);
% x0 = [id_l; id_q; id_e]; x = x0; y0 = zeros(m,1); y = y0; z0 = x0; z = z0;
% tau = 1; kappa = 1; theta0 = 1; theta = theta0;
% % Compute the auxiliary parameters which completely define (HSD)
% b_bar = (b*tau - A*x)/theta; 
% c_bar = (c*tau - A'*y - z)/theta; 
% g_bar = (c'*x - b'*y + kappa)/theta; 
% alpha_bar = (x'*z + tau*kappa)/theta;
% 
% % Formualte part of the coefficient matrix for the predictor search direction
% % Note that the system for the predictor search direction can be compactly written as
% % G_bar * dx_bar_p = R_bar_p, where G_bar = [G1; G2; G3] with
% G1 = [-A, sparse(m,m), sparse(m,Nt), b, sparse(m,1), -b_bar; sparse(Nt,Nt), A', speye(Nt), -c, sparse(Nt,1), c_bar; c', -b', zeros(1,Nt), 0, 1, -g_bar; -c_bar', b_bar', sparse(1,Nt), g_bar, 0, 0];
% % and x_bar is the concatenated vector of all variables
% x_bar = [x; y; z; tau; kappa; theta];
% 
% % Set the default relative accuracy
% rel_eps = 1e-8;
% max_iter_count = 100;
% display(['Algorithm begins with size(A) = [' num2str(size(A,1)) ',' num2str(size(A,2)) '] and density(A) = ' num2str(nnz(A)/(size(A,1)*size(A,2)))]);
% for it_cout =1:max_iter_count
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Form G_bar, R_bar_p and solve for the predictor search direction
%     G2 = [sparse(1,2*Nt+m), kappa, tau, 0]; G3 = [speye(Nt), sparse(Nt, m), theta*sparse(H_dual(z, dimension_info)), sparse(Nt,3)];
%     G_bar = [G1; G2; G3]; 
%     R_bar_p = [sparse(m+Nt+2,1); -tau*kappa; -x];
%     pred_dir = G_bar\R_bar_p; %linsolve(G_bar,R_bar_p);
%     dx_p = pred_dir(1:Nt); dy_p = pred_dir(Nt+1:Nt+m); dz_p = pred_dir(Nt+m+1: 2*Nt+m); dtau_p = pred_dir(2*Nt+m+1); dkappa_p = pred_dir(2*Nt+m+2); dtheta_p = pred_dir(2*Nt+m+3);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Compute alpha_p (approximately)
%     alpha_p = find_alpha_max(x_bar, pred_dir, dimension_info);
%     % Set sigma (parameter for the system for the combined search direction)
%     sigma = (1-alpha_p)^3;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Solve for the combined search direction
%     R_bar = [sparse(Nt+m+2, 1); -tau*kappa + sigma*theta - dtau_p*dkappa_p; -x - sigma*theta*g_dual(z, dimension_info)];
%     comb_dir = G_bar\R_bar; % linsolve(G_bar,R_bar);
%     dx = comb_dir(1:Nt); dy = comb_dir(Nt+1:Nt+m); dz = comb_dir(Nt+m+1:2*Nt+m); dtau = comb_dir(2*Nt+m+1); dkappa = comb_dir(2*Nt+m+2); dtheta = comb_dir(2*Nt+m+3);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Approximately find the step-length along the combined search direction and update the current iterate
%     alpha = 0.98*find_alpha_max(x_bar, comb_dir, dimension_info);
%     % x = x+alpha*dx; y = y + alpha*dy; z = z + alpha*dz; tau = tau + alpha*dtau; kappa = kappa + alpha*dkappa; theta = theta+alpha*dtheta;
%     x_bar = x_bar + alpha*comb_dir;
%     display(['theta=' num2str(theta,5) ',  sigma=' num2str(sigma,5) ',  dtheta=' num2str(dtheta,5) ',  alpha=' num2str(alpha,5) ', tau=' num2str(tau) ', kappa=' num2str(kappa)]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Termination conditions (from Section 5.4. in http://web.stanford.edu/~yyye/nonsymmhsdimp.pdf)
%     % Extract x, y, z, tau, theta, kappa
%     x = x_bar(1:Nt); y = x_bar(Nt+1:Nt+m); z = x_bar(Nt+m+1:2*Nt+m);  
%     tau = x_bar(2*Nt+m+1); kappa = x_bar(2*Nt+m+2); theta = x_bar(2*Nt+m+3);
%     % Check the 7 inequalities P, D, G, A, T, K, M
%     bool_P = norm(A*x-tau*b,Inf) <= rel_eps*max(1,norm([A,b],Inf)); 
%     bool_D = norm(A'*y+z-c*tau,Inf) <= rel_eps*max(1,norm([A',eye(Nt),-c], Inf)); 
%     bool_G = abs(-c'*x+b'*y-kappa) <= rel_eps*max(1, norm([-c',b',1],Inf)); 
%     bool_A = abs(c'*x/tau - b'*y/tau) <= rel_eps*(1+abs(b'*y/tau)); 
%     bool_T = tau <= rel_eps*(1e-2)*max(1,kappa); 
%     bool_K = tau <= rel_eps*(1e-2)*min(1,kappa); 
%     bool_M = theta <= rel_eps*(1e-2)*theta0;
%     if bool_P && bool_D && bool_A
%         x_final = x/tau; y_final = y/tau; z_final = z/tau;
%         display('The problem is primal and dual feasible.');
%         display(['The approximate primal and dual optimal objectives are ' num2str(c'*x/tau) ' and ' num2str(b'*y/tau)]);
%         % Generate x_return and z_return, cell arrays that is consistent with the input
%         idx_l = 1; idx_q = Nl+1; idx_e = Nl+sum(Nq)+1;
%         for k = 1:size(blk,1)
%             if blk{k,1} == 'l'
%                 x_return{k} = x_final(idx_l:idx_l+blk{k,2}-1); 
%                 z_return{k} = z_final(idx_l:idx_l+blk{k,2}-1); 
%                 idx_l = idx_l + blk{k,2};
%             elseif blk{k,1} == 'q'
%                 x_return{k} = x_final(idx_q:idx_q + sum(blk{k,2})-1);
%                 z_return{k} = z_final(idx_q:idx_q + sum(blk{k,2})-1); 
%                 idx_q = idx_q + sum(blk{k,2});
%             elseif blk{k,1} == 'e'
%                 x_return{k} = x_final(idx_e:idx_e + 3*blk{k,2}-1);
%                 z_return{k} = z_final(idx_e:idx_e + 3*blk{k,2}-1); 
%                 idx_e = idx_e + 3*blk{k,2};
%             end 
%         end
%         break;
%     elseif bool_P && bool_D && bool_G && bool_T
%         if c'*x < 0
%             display('The problem is dual infeasible and certificate of infeasibility has been found.');
%         elseif b'*y > 0
%             display('The problem is primal infeasible and a certificate of infeasibility has been found.');
%         end
%         break;
%     elseif bool_K && bool_M
%         display('The problem seems ill-posed!');
%         break;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Check whether the program gets stuck halfway
%     if abs(dtheta*alpha)/theta <= 1e-5
%         display('No more progress possible! The value of theta cannot decrease anymore! Algorithm terminates now!');
%         res1 = norm(A*x/tau-b, Inf)/norm([A, b], Inf); res2 = norm(A'*y/tau + z/tau -c, Inf)/norm([A' eye(Nt) c], Inf); res3 = abs(b'*y/tau - c'*x/tau)/(1+abs(b'*y/tau));
%         display(['Algorithm terminated with theta = ' num2str(theta) ', complementarity_residual = ' num2str((x'*z+tau*kappa-alpha_bar*theta)/(alpha_bar*theta0))]);
%         display(['The relative residual norms (linear primal, linear dual, duality gap) are ' num2str(res1) ', ' num2str(res2) ' and ' num2str(res3)]);
%         break;
%     end
% end