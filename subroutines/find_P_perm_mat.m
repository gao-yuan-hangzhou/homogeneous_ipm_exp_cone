function [P_perm_mat, all_removable_indices] = find_P_perm_mat(dense_col_indices, dimension_info)
% Given the indices of dense columns and dimension_info, 
% Find the appropriate permutation matrix such that
% A_hat*P = [A_hat_sp, A_hat_dc]
% where A_hat_dc contains all dense columns and some additional (possibly sparse) collumns 
% coming from the same blocks

Nl = dimension_info.l;
Nq = dimension_info.q;
Ne = dimension_info.e;
m = dimension_info.m;

n = Nl+sum(Nq)+3*Ne;

% Find all removable columns: flag them first
all_removable_columns_bool = zeros(n,1);
P_perm_mat = speye(n,n);

% If n is small, no need to find the dense columns
if n<1000
    all_removable_indices = [];
    return;
end

for k = 1:length(dense_col_indices)
    curr_dc_idx = dense_col_indices(k);
    if curr_dc_idx <= Nl
        % if dc_idx is a linear block, just move the current column
        all_removable_columns_bool(curr_dc_idx) = 1;
    elseif curr_dc_idx > Nl && curr_dc_idx <= Nl+sum(Nq)
        % For second-order cone part, move the whole block
        for idx_soc = 1:length(Nq)
            if (Nl+sum(Nq(1:idx_soc-1)) <= curr_dc_idx) && (curr_dc_idx <= Nl+sum(Nq(1:idx_soc)))
                curr_soc_indices = (Nl+sum(Nq(1:idx_soc-1)):Nl+sum(Nq(1:idx_soc)));
                break;
            end
        end
        all_removable_columns_bool(curr_soc_indices) = 1;
    else
        % For exp cone part, move all 3 columns
        curr_exp_idx = ceil((curr_dc_idx-(Nl+sum(Nq)))/3);
        all_removable_columns_bool(Nl+sum(Nq)+3*curr_exp_idx-2:Nl+sum(Nq)+3*curr_exp_idx) = 1;
    end
end

all_removable_indices = find(all_removable_columns_bool);
num_removable_columns = length(all_removable_indices);

% Begin constructing P_perm_mat
for k = 1:num_removable_columns
    removable_idx = all_removable_indices(k);
    new_idx = n - num_removable_columns + k;
    temp_col = P_perm_mat(:, removable_idx);
    P_perm_mat(:, removable_idx) = P_perm_mat(:, new_idx); 
    P_perm_mat(:, new_idx) = temp_col;
end

end

