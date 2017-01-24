% This is A_hat = A_hat_sp + A_hat_dc

% load('5th_iteration_data.mat');

% obtain the dense column indices

density_threshold = 0.1;

density = @(X) nnz(X)/numel(X);

dense_col_bool = sparse(dim_x, 1);

% A = Asp + Adc
A_hat_dc = sparse(m+2, dim_x);
A_hat_sp = A_hat;
for col_idx = 1:dim_x
    curr_col = A_hat(:, col_idx);
    if (density(curr_col) >= density_threshold)
        % disp([num2str(density(curr_col)), ', ', num2str(col_idx)]);
        dense_col_bool(col_idx) = 1;
        A_hat_dc(:, col_idx) = A_hat(:, col_idx);
        A_hat_sp(:, col_idx) = 0;
    end
end

dense_col_indices = find(dense_col_bool);

num_dc = length(dense_col_indices);
dc_A_hat = A_hat(:, dense_col_indices); % the matrix of all dense columns of A_hat

% construct the transformation matrix such that A_hat_dc = dc_A_hat*S_trans_mat
S_trans_mat = sparse(num_dc,dim_x);
for k = 1:num_dc
    dc_idx = dense_col_indices(k);
    S_trans_mat(k, dc_idx) = 1;
end

% Msp = A_hat_sp * H * A_hat_sp';
% Msp_pert = Msp+1e-8*norm(Msp, Inf)*speye(m+2);
% Lchol = chol(Msp_pert);

% Find UB, DB such that B_hat = UB*DB*UB'
U1 = [-b; 0; -g_bar/2]; U2 = [b_bar; g_bar/2; 0]; V1 = [zeros(m,1); 1; 0]; V2 = [zeros(m+1,1); 1];
UB = [U1, U2, V1, V2];
DB = [zeros(2), eye(2); -eye(2), zeros(2)];

% Now we consider A_hat = A_hat_sp + A_hat_dc 
% where A_hat_sp is the matrix of setting all dense columns to zero
% and A_hat_dc = dc_A_hat * S_trans_mat
% Let M = Msp + U*D*U' where
% Msp = A_hat_sp*H*A_hat_sp + A_hat_sp*H*A_hat_dc' + A_hat_dc*H*A_hat_sp' + diag([sparse(m,1); kappa/tau; 0])
% and U = [UB, dc_A_hat], D = blkdiag(DB, H_tilde)
% where H_tilde = S_trans_mat*H*S_trans_mat'

Msp = A_hat_sp*H*A_hat_sp' + A_hat_sp*H*A_hat_dc' + A_hat_dc*H*A_hat_sp' + diag([sparse(m,1); kappa/tau; 0]);
Msp_pert = Msp+1e-15*norm(Msp, Inf)*speye(m+2);
H_tilde = S_trans_mat*H*S_trans_mat';
U = [UB, dc_A_hat];
D = blkdiag(DB, H_tilde);

% By definition, 
% M = A_hat*H*A_hat' + B_hat + diag([sparse(m,1); kappa/tau; 0]);