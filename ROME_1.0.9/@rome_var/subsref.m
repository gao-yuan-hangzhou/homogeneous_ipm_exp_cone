function out = subsref(A,S)

% ROME_VAR\SUBSREF Implements subscripted referencing for rome_var objects
%
%   C1 = A(m)
%   C2 = A(m, n)
%
%
% Modification History: 
% 1. Joel 

% S

% default behavior if not using parentheses to subscript
if(~strcmp(S(1).type, '()'))
    out = builtin('subsref',A,S);
    return;
end

% if first subscript is parantheses,
S_index = S(1);

% number of subscripts
num_subscripts = length(S_index.subs);

% Case 1: Single index
if(num_subscripts == 1)

    % special case: colon is used (i.e for vectorization)
    if(strcmp(S_index.subs{1}, ':'))
        out = A;    % copy everything from A
        out.Size = [A.TotalSize, 1];
    else
        % standard single index case
        ind = S_index.subs{1};
        if(isvector(ind) && isvector(A))
            is_row_vec = (size(A, 1) == 1);
            if(is_row_vec)
                out = rome_var(1, length(ind));
            else
                out = rome_var(length(ind), 1);
            end
        else
            out = rome_var(size(ind));
        end
        out.BiAffineMap = A.BiAffineMap(ind(:), :); % select rows
        out.NumUnmappedVars     = A.NumUnmappedVars;
        out.NumUnmappedRandVars = A.NumUnmappedRandVars;
        out.NumMappedRandVars   = A.NumMappedRandVars;
        out.Cone = A.Cone;                          % preserve the continuity and cone constraints
        out.Continuity = A.Continuity;
        
        % special case: contiguous indexing - retain diagonal multiplier
        if(all(diff(ind(:)) == 1))
            % if A is diagonal, there may be some unnecssary
            % elements that I can (and should) remove
            if(A.IsDiag)
                num_pre_zero    = min(ind) - 1;
                out.BiAffineMap = [out.BiAffineMap(:, 1), out.BiAffineMap(:, 1 + (min(ind):max(ind)))];   % only select non-zero rows
                if(A.IsCertain)
                    out.NumUnmappedVars = out.NumUnmappedVars + num_pre_zero;                    
                elseif(A.IsRand)
                    out.NumUnmappedRandVars = out.NumUnmappedRandVars + num_pre_zero;
                    out.NumMappedRandVars = max(ind) - min(ind) + 1;
                else
                    error('rome_var:subsref:NoLDR', 'Do not allow Diagonal Optimization for LDR structures');
                end
            end
            
            % output has same diagonal multiplier as input
            out.DiagMult = A.DiagMult;
        end
    end

% Case 2: When we have multiple indices
else
    % check for illogical indices and return empty matrix
    if(any(cellfun('isempty', S_index.subs)))
        out = rome_empty_var;
        return;
    end

    % if all good, expand and vectorize indices
    S_index.subs = arrayfun(@expand_indices, S_index.subs, size(A), 'UniformOutput', false);
    dimension_array = cellfun('length', S_index.subs);

    % allocate output
    out = rome_var(dimension_array);

    % make a list of indices
    ind_list = cell(1, length(dimension_array));
    inner_reps_array = cumprod(dimension_array);
    inner_reps_array = [1, inner_reps_array];
    outer_reps_array = out.TotalSize ./ inner_reps_array(2:end);

    % iterate and make list
    for ii = 1:num_subscripts
        cur_subscripts = S_index.subs{ii};
        cur_ind_list = cur_subscripts(:, ones(1, inner_reps_array(ii))).';
        cur_ind_list = cur_ind_list(:);
        ind_list{ii} = repmat(cur_ind_list, outer_reps_array(ii), 1);
    end

    % get the index form of subscripts
    ind = sub2ind(size(A), ind_list{:});
    out.BiAffineMap = A.BiAffineMap(ind, :);
    out.NumUnmappedVars = A.NumUnmappedVars;
    out.NumUnmappedRandVars = A.NumUnmappedRandVars;
    out.NumMappedRandVars = A.NumMappedRandVars;
    out.Cone = A.Cone;                          % preserve the continuity and cone constraints
    out.Continuity = A.Continuity;

    % check contiguity
    if(all(diff(ind) == 1))
        out.DiagMult = A.DiagMult;
    end
    
    % squeeze out trailing singleton dimensions
    sz_out = size(out);
    for ii = length(sz_out):-1:3   %(stop at 3 - always at least a matrix)
        if(sz_out(ii) ~= 1)
            break;
        else
            out.Size(ii) = [];  % remove it
        end
    end
end

% if there is more than one index, we now call the built-in function
if(length(S) > 1)
    out = builtin('subsref',out,S(2:end));
end


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
