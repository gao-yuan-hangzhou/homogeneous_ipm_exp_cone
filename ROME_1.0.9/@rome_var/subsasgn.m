function A = subsasgn(A, S, B)

% ROME_VAR\SUBSASGN Implements subscripted assignment for rome_var objects
%
%
% Modification History: 
% 1. Joel 
% dbstack
% disp('SUBSASGN');

% default behavior if not using parentheses to subscript
if(~strcmp(S(1).type, '()'))
    A = builtin('subsasgn', A, S, B);
    return;
elseif(length(S) == 2)
    % doing a direct assignment to properties of the object (e.g. mean/cov)
    dummyvar = A.subsref(S(1));
    builtin('subsasgn', dummyvar, S(2:end), B);
    return;
end

% only want to override subscripted assignment on the primary object
% if first subscript is parantheses,
S_index = S(1);

% number of subscripts
num_subscripts = length(S_index.subs);

% if B is non-empty numeric - convert to rome-var
if(~isempty(B) && isnumeric(B))
    % convert constant matrices into rome_vars
    cur_val = B(:);       % get values of current input
    B = rome_var(size(B));
    B.BiAffineMap(:, 1) = cur_val;
    B.NumUnmappedVars = A.NumUnmappedVars;          % mimic A's unmapped entries
    B.NumUnmappedRandVars = A.NumUnmappedRandVars;
end

% Case 1: Single index
if(num_subscripts == 1)
    % special case: colon is used (assignment / deletion to all)
    if(strcmp(S_index.subs{1}, ':'))
        if(isempty(B))
            % define insertion index
            insert_ind = 1:A.TotalSize;
        else                
            % if B is non-empty, peform assignment
            % check on size
            if(B.TotalSize ~= A.TotalSize)
                error('rome_var:subsasgn:WrongNumberOfInputs', ...
                    'In an assignment  A(:) = B, the number of elements in A and B must be the same.');
            end

            % if B is scalar - expand it
            if(isscalar(B))
                B.BiAffineMap = repmat(B.BiAffineMap, A.TotalSize, 1);
                B.Size = [A.TotalSize, 1];
            end

            % overwrite A, retaining original size
            A = B; 
            return;
        end
    else
        % case of single index
        insert_ind = S_index.subs{:};
        
        % check for empty insert_inds
        if(isempty(insert_ind))
            return;
        end
        
        % don't allow indexing for matrices and N-D arrays
        if(~isvector(A) && max(insert_ind) > A.TotalSize)
            error('rome_var:subsasgn:CannotResizeArray', ...
                'In an assignment  A(I) = B, a matrix A cannot be resized.');
        end
        
        % compute new size of A
        num_expand = pos(max(insert_ind) - A.TotalSize);
        if(size(A, 1) == 1)
            A = cat(2, A, spalloc(1, num_expand, 0));   % if row vector
        else
            A = cat(1, A, spalloc(num_expand, 1, 0));   % if col vector
        end
    end
else
    % Case 2: N-dimensional subscripting
    % check for illogical indices and return the unchanged A
    if(any(cellfun('isempty', S_index.subs)))
        return;
    end

    % if all good, expand and vectorize indices
    S_index.subs = arrayfun(@expand_indices, S_index.subs, size(A), 'UniformOutput', false);
    dimension_array = cellfun('length', S_index.subs);

    % make a list of indices
    ind_list = cell(1, length(dimension_array));
    inner_reps_array = [1, cumprod(dimension_array)];
    outer_reps_array = inner_reps_array(end) ./ inner_reps_array(2:end);

    % iterate and make list
    for ii = 1:num_subscripts
        cur_subscripts = S_index.subs{ii};
        cur_ind_list = cur_subscripts(:, ones(1, inner_reps_array(ii))).';
        cur_ind_list = cur_ind_list(:);
        ind_list{ii} = repmat(cur_ind_list, outer_reps_array(ii), 1);
    end
    
    % Expand A
    max_ind = cellfun(@max, ind_list);
    expand_sz = max(size(A), max_ind) - size(A);
    for ii = 1:numel(expand_sz)        
        if(expand_sz == 0)
            continue;
        end
        
        cur_sz = A.size;   % notice that size (A) changes on each iteration
        cur_sz(ii) = expand_sz(ii);
        A = cat(ii, A, zeros(cur_sz));
    end

    % get the index form of subscripts
    insert_ind = sub2ind(A.Size, ind_list{:});
end

% Handle Deletion
if(isempty(B))   
    % Index of non-colons
    numeric_ind = cellfun(@isnumeric, S(1).subs);
    
    % if too many numerics,
    if(sum(numeric_ind) > 1)
        error('rome_var:subsasgn:TooManyNumeric', 'A null assignment can have only one non-colon index.');
    % if all are colons
    elseif(~any(numeric_ind))
        A = [];
    else
        % if single non-colon, standard Deletion
        A.BiAffineMap(insert_ind, :) = [];
        A.DiagMult = NaN;

        % Single Subcript case
        if(num_subscripts == 1)
            if(size(A, 2) == 1)
                % column vectors
                A.Size = [size(A, 1) - numel(insert_ind), 1];
            else
                % matrices, N-D arrays and row-vectors handled similarly
                A.Size = [1, A.TotalSize - numel(insert_ind)];
            end

        % Multiple Subscript Case
        else
            oldsz_A = A.Size;
            oldsz_A(numeric_ind) = oldsz_A(numeric_ind) - numel(S_index.subs{numeric_ind});
            A.Size = oldsz_A;
        end
    end
    
    % Finish Deletion, return
    return;
end

% if B is scalar - expand it
if(isscalar(B))
    B.BiAffineMap = repmat(B.BiAffineMap, numel(insert_ind), 1);
    B.Size = [numel(insert_ind), 1];
end

% PERFORM ASSIGNMENT
% -------------------
% define input arguments
first_map = A.BiAffineMap;
second_map = B.BiAffineMap;

% just assign.
num_pre_zero = A.NumUnmappedVars - B.NumUnmappedVars;
num_post_zero= A.NumMappedVars - B.NumMappedVars + num_pre_zero;
num_pre_zero_rand = A.NumUnmappedRandVars  - B.NumUnmappedRandVars;
num_post_zero_rand = A.NumMappedRandVars  - B.NumMappedRandVars + num_pre_zero_rand;

% optimization to reduce extraneous dependant vars
if(B.IsRand)
    num_pre_zero = -A.NumMappedVars;
    num_post_zero = 0;
end
if(B.IsCertain)
    num_pre_zero_rand = -A.NumMappedRandVars;
    num_post_zero_rand = 0;
end
if(A.IsRand)
    num_pre_zero = B.NumMappedVars;
    num_post_zero = 0;
    A.NumUnmappedVars = B.NumUnmappedVars;
end
if(A.IsCertain)    
    num_pre_zero_rand = B.NumMappedRandVars;
    num_post_zero_rand = 0;
    A.NumUnmappedRandVars = B.NumUnmappedRandVars;
end


% get output number of cols
inner_sz = (1 + A.NumMappedVars) + pos(num_pre_zero) + neg(num_post_zero);
out_num_cols = inner_sz .* ((1 + A.NumMappedRandVars) + pos(num_pre_zero_rand) + neg(num_post_zero_rand));

% for convenience, get number of rows (corresponds to number of vars)
num_rows = A.TotalSize;

% allocate memory
expanded_first_map = spalloc(num_rows, out_num_cols, nnz(first_map));
first_inner_ind = [1, pos(num_pre_zero) + 1 + (1:A.NumMappedVars)];     % inner indices

% make indices by summing
first_inner_ind_expand = outer_sum(first_inner_ind', 0:inner_sz:(inner_sz*(A.NumMappedRandVars-1)));
first_inner_ind_expand = first_inner_ind_expand(:)';

% make the complete index vector
first_ind_vec = [first_inner_ind, inner_sz * (pos(num_pre_zero_rand) + 1) + first_inner_ind_expand(:)'];

% if first map is empty, don't assign
if(~isempty(first_map))
    % assign
    expanded_first_map(:, first_ind_vec) = first_map;
end

% allocate memory
expanded_second_map = spalloc(size(second_map, 1), out_num_cols, nnz(second_map));
second_inner_ind = [1, neg(num_pre_zero) + 1 + (1:B.NumMappedVars)];     % inner indices

% make indices by summing
second_inner_ind_expand = outer_sum(second_inner_ind', 0:inner_sz:(inner_sz*(B.NumMappedRandVars-1)));
second_inner_ind_expand = second_inner_ind_expand(:)';

% make the complete index vector
second_ind_vec = [second_inner_ind, inner_sz * (neg(num_pre_zero_rand) + 1) + second_inner_ind_expand(:)' ];

% assign
expanded_second_map(:, second_ind_vec) = second_map;

% make subscripted assigment
expanded_first_map(insert_ind, :) = expanded_second_map;

% perform assignment
A.BiAffineMap = expanded_first_map;
A.DiagMult = NaN;     % assignment usually destroys diagonality.

% Assign number of unmapped variables
A.NumUnmappedVars = min(A.NumUnmappedVars, B.NumUnmappedVars);
A.NumUnmappedRandVars = min(A.NumUnmappedRandVars,B.NumUnmappedRandVars);

% define number of mapped variables
A.NumMappedRandVars  = (size(A.BiAffineMap, 2) ./ inner_sz) - 1;

%
% EXPAND_INDICES Helper function to expand and vectorize indices
%
function ind = expand_indices(ind, max_sz)
if(numel(ind) == 1 && strcmp(ind, ':'))
    ind = (1:max_sz)';
else
    ind = ind{:};
    ind = ind(:);
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
