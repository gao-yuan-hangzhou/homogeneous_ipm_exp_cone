function H_overall = H_exp_tilde_dual_sparse_diagonal(input_vec)
% Compute the Hessian of the dual exp barrier (from the second pair)
% where N is the number of dual exp cones
Ne = length(input_vec)/3;
if floor(Ne)~=Ne
    display('Error: in g_exp_tilde_dual, the input vector has dimension not equal to 3*Ne');
end

% For large Ne, we do the following: 
% partition it into smaller (dense) blocks, say, 9000-by-9000, i.e., 3000 exp blocks a time
Ne_per_block = 3000;
num_dense_blocks = ceil(Ne/Ne_per_block);
last_dense_block_size = Ne - Ne_per_block*(num_dense_blocks-1);

H_overall = sparse([]);
for blk_idx = 1:num_dense_blocks
    H_curr = zeros(3*Ne_per_block, 3*Ne_per_block);
    kmax = Ne_per_block;
    if blk_idx == num_dense_blocks
        kmax = last_dense_block_size;
    end
    for k = 1:kmax
        curr_exp_idx = Ne_per_block*(blk_idx-1)+k;
        u = input_vec(3*curr_exp_idx-2); v = input_vec(3*curr_exp_idx-1); w = input_vec(3*curr_exp_idx);
        r = log(-v/u); t = u-w+u*r;
        H_curr(3*k-2:3*k, 3*k-2:3*k) = [      1/(u*t) + r^2/t^2 + 1/u^2,             u*r/(v*t^2) - 1/(v*t),         -r/t^2;
                                        (u*log(-v/u))/(v*t^2) - 1/(v*t), 1/v^2 + u/(v^2*t) + u^2/(v^2*t^2),     -u/(v*t^2);
                                                                 -r/t^2,                        -u/(v*t^2),          1/t^2];
    end
    H_overall = blkdiag(H_overall, H_curr);
end
    
% for k=1:Ne
%     u = input_vec(3*k-2); v = input_vec(3*k-1); w = input_vec(3*k); 
%     r = log(-v/u); t = u-w+u*r;
%     % Evaluate the Hessian of f_tilde_dual,
%     % The (second) barrier function for the dual of the exponential cone
%      H_overall(3*k-2:3*k, 3*k-2:3*k) = [      1/(u*t) + r^2/t^2 + 1/u^2,             u*r/(v*t^2) - 1/(v*t),         -r/t^2;
%                                 (u*log(-v/u))/(v*t^2) - 1/(v*t), 1/v^2 + u/(v^2*t) + u^2/(v^2*t^2),     -u/(v*t^2);
%                                                          -r/t^2,                        -u/(v*t^2),          1/t^2];
%     
end

