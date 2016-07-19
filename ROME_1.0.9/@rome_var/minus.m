function out_var_obj = minus(A, B)

% ROME_VAR\MINUS Implements matrix subtraction in ROME
%
%   C = A - B 
%
%   At least one of A, B must be a rome_var .
%
% Modification History: 
% 1. Joel 

% operation on an empty matrix will yield the other
if(isempty(A))
    out_var_obj = -B;
    return;
elseif(isempty(B))
    out_var_obj = A;
    return;
end


% Throw an error if the sizes do not match.
if(~isscalar(A) && ~isscalar(B) && any(size(A) ~= size(B)))
    error('rome_var:minus:InvalidArg', 'Input matrices must have the same size.');
end

% Allocate output variable
if(isscalar(A))
    out_var_obj = rome_var(size(B));
else
    out_var_obj = rome_var(size(A));
end


% Case 1: A is a (constant) numeric matrix and B is a rome_var
if(isnumeric(A))
    if(isscalar(B))
        out_var_obj.BiAffineMap = [A(:) - B.BiAffineMap(:, 1), ...
                                    repmat(-B.BiAffineMap(:, 2:end), numel(A), 1)];
    else
        out_var_obj.BiAffineMap = [A(:) - B.BiAffineMap(:, 1), -B.BiAffineMap(:, 2:end)];
    end
    out_var_obj.DiagMult = -B.DiagMult;
    out_var_obj.NumUnmappedVars = B.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = B.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars  = B.NumMappedRandVars;
    
% Case 2: A is a rome_var and B is a (constant) double matrix
elseif(isnumeric(B))
    if(isscalar(A))
        out_var_obj.BiAffineMap = [A.BiAffineMap(:, 1) - B(:), ...
                                    repmat(A.BiAffineMap(:, 2:end), numel(B), 1)];
    else
        out_var_obj.BiAffineMap = [A.BiAffineMap(:, 1) - B(:), A.BiAffineMap(:, 2:end)];
    end
    out_var_obj.DiagMult = A.DiagMult;
    out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars  = A.NumMappedRandVars;


% Case 3: Both are rome_vars
else
    % define input arguments
    first_map = A.BiAffineMap;
    second_map = B.BiAffineMap;
    
    % expand scalars
    if(isscalar(A))
        first_map = repmat(first_map, prod(size(B)), 1);
    elseif(isscalar(B))
        second_map = repmat(second_map, prod(size(A)), 1);
    end

    % Assign number of unmapped variables
    out_var_obj.NumUnmappedVars = min(A.NumUnmappedVars, B.NumUnmappedVars);
    out_var_obj.NumUnmappedRandVars = min(A.NumUnmappedRandVars, B.NumUnmappedRandVars);
   
    % calculate number of required zero-pads
    num_pre_zero = A.NumUnmappedVars - B.NumUnmappedVars;
    num_post_zero= A.NumMappedVars - B.NumMappedVars + num_pre_zero;
    num_pre_zero_rand = A.NumUnmappedRandVars  - B.NumUnmappedRandVars;
    num_post_zero_rand = A.NumMappedRandVars  - B.NumMappedRandVars + num_pre_zero_rand;
    
    % optimization to reduce extraneous dependant vars
    if(B.IsRand)
        num_pre_zero = -A.NumMappedVars;
        num_post_zero = 0;
        out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
    end
    if(B.IsCertain)
        num_pre_zero_rand = -A.NumMappedRandVars;
        num_post_zero_rand = 0;
        out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
    end    
    
    if(A.IsRand)
        num_pre_zero = B.NumMappedVars;
        num_post_zero = 0;
        out_var_obj.NumUnmappedVars = B.NumUnmappedVars;
    end
    if(A.IsCertain)
        num_pre_zero_rand = B.NumMappedRandVars;
        num_post_zero_rand = 0;
        out_var_obj.NumUnmappedRandVars = B.NumUnmappedRandVars;
    end
    
    % get output number of cols
    inner_sz = (1 + A.NumMappedVars) + pos(num_pre_zero) + neg(num_post_zero);
    out_num_cols = inner_sz .* ((1 + A.NumMappedRandVars) + pos(num_pre_zero_rand) + neg(num_post_zero_rand));

    % for convenience, get number of rows (corresponds to number of vars)
    num_rows = out_var_obj.TotalSize;
    
    % allocate memory
    expanded_first_map = spalloc(num_rows, out_num_cols, nnz(first_map));    
    first_inner_ind = [1, pos(num_pre_zero) + 1 + (1:A.NumMappedVars)];     % inner indices
    
    % make indices by summing
    first_inner_ind_expand = outer_sum(first_inner_ind', 0:inner_sz:(inner_sz*(A.NumMappedRandVars-1)));
     
    % make the complete index vector
    first_ind_vec = [first_inner_ind, inner_sz * (pos(num_pre_zero_rand) + 1) + first_inner_ind_expand(:)'];
                 
    % assign
    expanded_first_map(:, first_ind_vec) = first_map;
                 
    % allocate memory
    expanded_second_map = spalloc(num_rows, out_num_cols, nnz(second_map));
    second_inner_ind = [1, neg(num_pre_zero) + 1 + (1:B.NumMappedVars)];     % inner indices
    
    % make indices by summing
    second_inner_ind_expand = outer_sum(second_inner_ind', 0:inner_sz:(inner_sz*(B.NumMappedRandVars-1)));
        
    % make the complete index vector
    second_ind_vec = [second_inner_ind, inner_sz * (neg(num_pre_zero_rand) + 1) + second_inner_ind_expand(:)' ];
    
    % assign
    expanded_second_map(:, second_ind_vec) = second_map;
    
    % perform zero-padded subtraction
    out_var_obj.BiAffineMap = expanded_first_map - expanded_second_map;

    % define number of unmapped variables
    out_var_obj.NumMappedRandVars  = (size(out_var_obj.BiAffineMap, 2) ./ inner_sz) - 1;
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
