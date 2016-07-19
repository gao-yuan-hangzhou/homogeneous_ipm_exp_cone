function out_var_obj = sum(A, dim)

% ROME_VAR\SUM Sums up collection of rome_vars
%
%   C = sum(A, dim)
%
%   C = sum(A)
%
% Modification History: 
% 1. Joel 

% catch pathological case
if(isempty(A))
    out_var_obj = [];
    return;
end

% get original size
orig_size = size(A);

% handle default for second argument
if(nargin < 2) 
    dim = find(size(A) > 1);
    if(~isempty(dim))
        dim = dim(1);
    else
        dim = 1;
    end
elseif(dim < 1 || ~isscalar(dim) || floor(dim) ~= dim)
    error('rome_var:sum:InvalidArg', ...
        'Summing dimension must be positive scalar integer within indexing range');
elseif(dim > length(orig_size))
    % summing over a larger dimesion returns yourself
    out_var_obj = A;
    return;
end

% create new size vector, squeeze out trailing singleton dimensions
new_size = [orig_size(1:(dim-1)), 1, orig_size((dim+1):end)];
if(new_size(end) == 1)
    new_size = new_size(1:end-1);
end

% allocate memory for new object
out_var_obj = rome_var(new_size);
out_var_obj.NumUnmappedVars = A.NumUnmappedVars; % no change to unmapped variables
out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars; % no change to unmapped variables
out_var_obj.NumMappedRandVars = A.NumMappedRandVars; % no change to unmapped variables

% create indexing array
num_to_sum = orig_size(dim);
num_input_rows = A.TotalSize;
num_output_rows = prod(new_size);
inner_dim = prod(orig_size(1:(dim-1)));

% ref for matrix permutation:
% http://home.online.no/~pjacklam/matlab/doc/mtt/doc/archive/2000-05-05/mtt.pdf.gz
ind = reshape((1:num_input_rows)', inner_dim, num_input_rows/ inner_dim);
col_ind = reshape(ind, inner_dim, 1, num_to_sum, num_input_rows / num_to_sum / inner_dim);
col_ind = permute(col_ind, [1, 4, 3, 2]);
col_ind = reshape(col_ind, num_output_rows, num_to_sum); 

% construct the indices
row_ind = repmat((1:num_output_rows)', 1, num_to_sum);

% make the summation matrix
sum_matrix = sparse(row_ind(:), col_ind(:), 1, num_output_rows, num_input_rows, num_input_rows);

% multiply in
out_var_obj.BiAffineMap = sum_matrix*A.BiAffineMap;

% % make the summation matrix
% sum_matrix = spalloc(num_output_rows, num_input_rows, num_input_rows);
% sum_matrix(sub2ind(size(sum_matrix), row_ind(:), col_ind(:))) = 1;
% 
% % multiply in
% out_var_obj.BiAffineMap = sum_matrix*A.BiAffineMap;



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
