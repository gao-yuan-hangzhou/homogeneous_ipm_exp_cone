% Comparison of the 3 algorithms
% ======================
% Defining the assignment problem min <C,X> s.t. Xe = a, X'e = b

addpath(genpath(pwd))
rng('default')
DIM = 100;
m = DIM; n = DIM; C = rand(m,n); a = ones(m,1); b = ones(n,1);
C = C/10;
total_number_of_iterations = 10000;

% ======================
% Formulate the assignment problem into standard form LP min c'x s.t. Ax=b
A = [kron(ones(1,m), speye(n)); kron(speye(n),ones(1,m))]; bb = [a;b];
% Truncate A and b so that A is full row rank
A = A(1:end-1,:); bb = bb(1:end-1); 
cc = reshape(C', [m*n,1]);

% ======================
% Solve the LP using Gurobi
if (true)
   clear model; 
   model.obj = cc; model.A = sparse(A); model.rhs = bb; model.sense = '='; 
   clear params;
   params.Presolve = 2; params.TimeLimit = 100;
   result = gurobi(model, params);  
   xg = result.x; Xg = reshape(result.x, [m,n]); optimal_cost = result.objval;
end
% ======================
% Solve the LP using Bregman ADMM (Wang & Banerjee)
if (true)
    tstart = clock;
    [x1, badmm_optimal_obj_val, badmm_iteration_count,badmm_res] ...
    = badmm_assignment_problem_solver(C, a, b, 0.1, total_number_of_iterations);
    ttime = etime(clock,tstart);
    fprintf('\n badmm_optimal_obj_val = %7.6f',badmm_optimal_obj_val); 
    fprintf('\n badmm_iteration_count = %2.0f',badmm_iteration_count); 
    fprintf('\n badmm_iteration_time  = %3.2f\n',ttime);
    hold off;
    semilogy(badmm_res); hold on;
end

% ======================
% Solve the LP using Semi-proximal ALM
tstart = clock;
[x2,y2,z2,res2] = semi_proximal_alm_lp_solver(cc, A, bb);
ttime = etime(clock,tstart);
fprintf('\n spalm_optimal_obj_val = %7.6f',cc'*x2); 
fprintf('\n spalm_iteration_count = %2.0f',length(res2));
fprintf('\n spalm_iteration_time  = %3.2f\n',ttime);
semilogy(res2,'r');
% ======================
% Solve the LP using ADMM (with He & Yuan's formulation)
if (n <= 100)
   tstart = clock;
   [x3,y3,z3,res3] = admm_reformulated(cc, A, bb);  
   ttime = etime(clock,tstart);
   fprintf('\n admm_optimal_obj_val  = %7.6f',cc'*x3); 
   fprintf('\n admm_iteration_count  = %2.0f\n',length(res3));
   fprintf('\n admm_iteration_time   = %3.2f\n',ttime);   
   semilogy(res3, 'g')
end
fprintf('\n')
