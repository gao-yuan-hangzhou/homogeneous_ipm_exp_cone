% This is A_hat*P = [A_hat_sp, A_hat_dc]

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

U1 = [-b; 0; -g_bar/2]; U2 = [b_bar; g_bar/2; 0];

V1 = [zeros(m,1); 1; 0]; V2 = [zeros(m+1,1); 1];

UB = [U1, U2, V1, V2];
DB = [zeros(2), eye(2); -eye(2), zeros(2)];

% By definition, the Schur complement matrix
% M = A_hat * H * A_hat' + B_hat + diag([sparse(m,1); kappa/tau; 0])
% where B_hat can be written as
% B_hat = UB * DB * UB'
% Assume A_hat has s dense columns with dense_col_indices = [i(1), ..., i(s)]
% To handle the dense columns in A_hat, let P be an appropriate permutation such that
% A_hat*P = [A_hat_sp, A_hat_dc], or equivalently, 
% A_hat = [A_hat_sp, A_hat_dc]*P'
% where A_hat_dc = A_hat(:, dense_col_indices) and 
% A_hat_sp contains the remaining columns of A_hat.
% Then A_hat * H * A_hat' = A_hat_sp * H_tilde_1 * A_hat_sp' + A_hat_dc * H_tilde_2 * A_hat_dc'
% where H_tilde = P'*H*P,
% H_tilde_1 = H_tilde(1:n-s, 1:n-s); H_tilde_2 = H_tilde(n-s+1:end, n_s+1:end)
% As such, we have 
% Msp = A_hat_sp * H_tilde_1 * A_hat_sp' + diag([sparse(m,1); kappa/tau; 0])
% U = [UB, A_hat_dc];
% D = blkdiag(DB, H_tilde_2);

% Construct the permutation: 
n = size(A_hat,2);
P_perm_mat = speye(n, n);
% A_hat_copy = A_hat;
for k = 1:num_dc
    dc_idx = dense_col_indices(k);
    new_idx = n - num_dc + k;
    temp_col = P_perm_mat(:, dc_idx);
    P_perm_mat(:, dc_idx) = P_perm_mat(:, new_idx); 
    P_perm_mat(:, new_idx) = temp_col;
end

% Construct A_hat_sp and A_hat_dc
A_hat_sp = A_hat * P_perm_mat(:, 1:n-num_dc);
A_hat_dc = A_hat * P_perm_mat(:, n-num_dc+1:end);

% Construct H_tilde, H_tilde_1, H_tilde_2
H_tilde = P_perm_mat'*H*P_perm_mat;
H_tilde_1 = P_perm_mat(:, 1:n-num_dc)'*H*P_perm_mat(:, 1:n-num_dc);
H_tilde_2 = P_perm_mat(:, n-num_dc+1:end)'*H*P_perm_mat(:, n-num_dc+1:end);

% Set Msp, U, D such that M = Msp + U*D*U';
Msp = A_hat_sp * H_tilde_1 * A_hat_sp' + diag([sparse(m,1); kappa/tau; 0]);
Msp_pert = Msp+1e-15*norm(Msp, Inf)*speye(m+2);
U = [UB, A_hat_dc];
D = blkdiag(DB, H_tilde_2);
% By definition:
% M = A_hat * H * A_hat' + diag([sparse(m,1); kappa/tau; 0]) + B_hat;
% Check the linear algebra:
% norm(M - (Msp + U*D*U'), Inf)
