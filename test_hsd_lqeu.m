addpath ./subroutines % All subroutines are in the folder "subroutines"

% We generate (feasible) instances of conic programs with 
% linear, second-order and exponential cone constraints
% to test hsd_lqe.m
% Choose whether to assure feasibility of the generated instance
ensure_feasibility = true;
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
clear blk; clear c_cell; clear A_cell; m = 18;
blk{1,1} = 'u'; blk{1,2} = 2; A_cell{1} = sprandn(m,sum(blk{1,2}),0.05); c_cell{1} = randn(blk{1,2},1);
blk{2,1} = 'l'; blk{2,2} = 800; A_cell{2} = sprandn(m, sum(blk{2,2}), 0.1); c_cell{2} = randn(sum(blk{2,2}),1);
blk{3,1} = 'e'; blk{3,2} = 3*ones(50,1); A_cell{3} = sprandn(m, sum(blk{3,2}), 0.15); c_cell{3} = randn(sum(blk{3,2}), 1);
blk{4,1} = 'q'; blk{4,2} = max(2,randi(50,5,1)); A_cell{4} = sprandn(m, sum(blk{4,2}), 0.12); c_cell{4} = randn(sum(blk{4,2}), 1);
%blk{5,1} = 'u'; blk{5,2} = 3; A_cell{5} = sprandn(m, blk{5,2}, 0.08); c_cell{5} = randn(blk{5,2}, 1);
%blk{6,1} = 'l'; blk{6,2} = 1; A_cell{6} = sprandn(m, blk{6,2}, 0.07); c_cell{6} = randn(blk{6,2}, 1);
%blk{7,1} = 'e'; blk{7,2} = 3*ones(10,1); A_cell{7} = sprandn(m, sum(blk{7,2}), 0.11); c_cell{7} = randn(sum(blk{7,2}),1);
b = randn(m,1);

% Choose whether to have a feasible instance
if ensure_feasibility
    clear A_cell c_cell;
    [A_cell, c_cell, b] = generate_random_feasible_instance(blk,m);
end

%save('blk_input.mat', 'blk', 'A_cell', 'c_cell', 'b');
load('blk_input.mat', 'blk', 'A_cell', 'c_cell', 'b');

%[obj_val, x_return,y_return,z_return, info] = hsd_lqeu_NT_Mehrotra(blk, A_cell, c_cell, b, 1e-8, 1000);
[obj_val, x_return,y_return,z_return, info] = hsd_lqeu(blk, A_cell, c_cell, b, 1e-8, 1000);

% Check certificate of dual infeasibility
if strcmp(info.solution_status, 'dual_infeasible') || strcmp(info.solution_status, 'primal_and_dual_infeasible')
    cTx = 0;
    for k = 1:size(blk,1)
        cTx = cTx + c_cell{k}' * info.x_certificate_dual_infeasible{k};
    end
    disp(['cert_dual_infeas = cTx = ' num2str(cTx) ' < 0']);
end

% Check certificate of primal infeasibility
if strcmp(info.solution_status, 'primal_infeasible') || strcmp(info.solution_status, 'primal_and_dual_infeasible')
    bTy = b'*info.y_certificate_primal_infeasible;
    disp(['cert_primal_infeas = bTy = ' num2str(bTy) ' > 0']);
end
% Check objective value
if strcmp(info.solution_status, 'optimal')
    bTy = b'*y_return; cTx = 0;
    for k = 1:size(blk,1)
        cTx = cTx + c_cell{k}' * x_return{k};
    end
    disp('Check returned optimal solutions in cell arrays.');
    display(['bTy = ' num2str(bTy)]);
    display(['cTx = ' num2str(cTx)]);
end