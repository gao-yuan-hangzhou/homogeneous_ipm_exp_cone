addpath ./subroutines
run ./SDPT3-4.0/startup.m;

% We generate (feasible) instances of conic programs with 
% LINEAR and SECOND-ORDER cone constraints
% to test hsd_lqe.m against sdpt3.m
% Choose whether to assure feasibility of the generated instance

assure_feasible_instance = true;

clear blk; clear c_cell; clear A_cell; 
% The total number of linear constraints
% The different types of blocks
m = 50;
blk{1,1} = 'l'; blk{1,2} = 100; A_cell{1} = sprandn(m,sum(blk{1,2}), 0.12); c_cell{1} = randn(sum(blk{1,2}),1);
blk{2,1} = 'u'; blk{2,2} = 1; A_cell{2} = sprandn(m,sum(blk{2,2}), 0.18);  c_cell{2} = zeros(sum(blk{2,2}),1);
blk{3,1} = 'q'; blk{3,2} = max(2,randi(50,5,1)); A_cell{3} = sprandn(m,sum(blk{3,2}), 0.15);  c_cell{3} = randn(sum(blk{3,2}),1);
blk{4,1} = 'u'; blk{4,2} = 1; A_cell{4} = sprandn(m,sum(blk{4,2}), 0.5); c_cell{4} = zeros(sum(blk{4,2}),1);
% The random RHS vector for Ax=b, where A = [A(1), ..., A(N)], x = [x(1); ...; x(N)]
b = randn(m,1);

if assure_feasible_instance
    clear A_cell c_cell b;
    [A_cell, c_cell, b] = generate_random_feasible_instance(blk,m);
end

% Solve the instance using SDPT3
[obj_sdpt3,x_sdpt3,y_sdpt3,z_sdpt3,info_sdpt3,runhist_sdpt3] = sdpt3(blk,A_cell,c_cell,b);

% Solve the instance using hsd_lqe
[obj_lqe, x_lqe, y_lqe, z_lqe, info_lqe] = hsd_lqeu(blk, A_cell, c_cell,b);

% Compare the (primal and dual) objective values
display(obj_sdpt3);
display(obj_lqe);