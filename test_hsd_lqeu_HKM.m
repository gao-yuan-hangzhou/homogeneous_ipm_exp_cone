addpath ./subroutines
ensure_feasibility = true;
clear blk; clear c_cell; clear A_cell; m = 18;
blk{1,1} = 'u'; blk{1,2} = 2; A_cell{1} = sprandn(m,sum(blk{1,2}),0.05); c_cell{1} = randn(blk{1,2},1);
blk{2,1} = 'l'; blk{2,2} = 800; A_cell{2} = sprandn(m, sum(blk{2,2}), 0.1); c_cell{2} = randn(sum(blk{2,2}),1);
blk{3,1} = 'e'; blk{3,2} = 3*ones(50,1); A_cell{3} = sprandn(m, sum(blk{3,2}), 0.15); c_cell{3} = randn(sum(blk{3,2}), 1);
blk{4,1} = 'q'; blk{4,2} = max(2,randi(50,5,1)); A_cell{4} = sprandn(m, sum(blk{4,2}), 0.12); c_cell{4} = randn(sum(blk{4,2}), 1);
if ensure_feasibility
    clear A_cell c_cell; [A_cell, c_cell, b] = generate_random_feasible_instance(blk,m);
end

save('blk_input.mat', 'blk', 'A_cell', 'c_cell', 'b');
load('blk_input.mat', 'blk', 'A_cell', 'c_cell', 'b');

hsd_lqeu_HKM(blk, A_cell, c_cell, b);