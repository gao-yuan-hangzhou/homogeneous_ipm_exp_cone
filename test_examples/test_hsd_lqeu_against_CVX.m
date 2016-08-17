addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']); addpath([fileparts(pwd), '/cvx']);
% cvx_setup;
% clear;

% Set the number of exponential cones, i.e. dimension of x = 3*Ne
Ne = 5000; Nl = 150; Nq = max(2,randi(50,10,1));

% Construct a feasible instance
blk{1,1} = 'e'; blk{1,2} = 3*ones(Ne,1); 
blk{2,1} = 'l'; blk{2,2} = Nl;
blk{3,1} = 'q'; blk{3,2} = Nq;
[A_cell, c_cell, b] = generate_random_feasible_instance(blk,1);

% Solve the problem using hsd_lqeu
[opt_sol, x_retun, y_return, z_return, info] = hsd_lqeu_Schur(blk, A_cell, c_cell, b);
disp(['dual optimal solution by hsd_lueq = ' num2str(opt_sol(1))]);

% Solve the problem using CVX
cvx_clear
cvx_tic
cvx_begin
variable xe(3*Ne)
variable xl(Nl)
variable xq(sum(Nq))
minimize(c_cell{1}'*xe + c_cell{2}'*xl+c_cell{3}'*xq)
subject to
    A_cell{1}*xe + A_cell{2}*xl + A_cell{3}*xq == b;
    xl >= 0
    for k = 1:Ne
        {xe(3*k-2), xe(3*k), xe(3*k-1)} <In> exponential
    end
    for k = 1:length(Nq)
        {xq(sum(Nq(1:k-1))+2:sum(Nq(1:k))), xq(sum(Nq(1:k-1))+1)} <In> lorentz(Nq(k)-1)
    end
cvx_end
cvx_toc