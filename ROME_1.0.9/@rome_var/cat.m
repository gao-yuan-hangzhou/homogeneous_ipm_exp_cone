function out_var_obj = cat(dim, varargin)

% ROME_VAR\CAT Implements N-D concatenation for rome_var
%
%   C = cat(3, A, B .. );
%
%   TEMP: Still in progress! Only works for very simple case
%
% Modification History: 
% 1. Joel 


% define number of rome_var objects here
num_objs = numel(varargin);

% Optimization: if nargin = 2 just return
if(num_objs == 0)
    error('rome_var:cat:MissingObject', 'Need to supply at least 1 rome_var object to concatenate');
elseif(num_objs == 1)
    out_var_obj = varargin{1};
    return;
end

% Note that nargin represents how many matrices we are concatenating
% parse arguments and extract size information
num_mapped_var_arr   = zeros(1, num_objs);
num_unmapped_var_arr = zeros(1, num_objs);
num_mapped_rand_var_arr   = zeros(1, num_objs);
num_unmapped_rand_var_arr = zeros(1, num_objs);
is_pure_rand_arr = zeros(1, num_objs);
is_pure_certain_arr = zeros(1, num_objs);
cat_dim_arr = zeros(1, num_objs);
sz_dims = [];

for ii = 1:num_objs
    cur_obj = varargin{ii};
    
    % ignore empty matrices
    if(isempty(cur_obj))
        is_pure_rand_arr(ii) = 1;
        is_pure_certain_arr(ii) = 1;
        continue;
    end
    
    % convert constant matrices into rome_vars
    if(isnumeric(cur_obj))
        cur_val = cur_obj(:);       % get values of current argument
        cur_obj = rome_var(size(cur_obj));
        cur_obj.BiAffineMap(:, 1) = cur_val;
        varargin{ii} = cur_obj;     % re-assign argument
    end
    
    % compute object size
    sz_obj = size(cur_obj);
    if(numel(sz_obj) < dim)    % expand trailing singleton dimensions
        sz_obj((numel(sz_obj) + 1):dim) = 1;
    end

    % for catting dimension, output size is sum of sizes
    cat_dim_arr(ii) = sz_obj(dim);
    
    % first case, assign sz obj
    if(isempty(sz_dims))
        sz_dims = sz_obj;
    else
        % check that the rest of the sizes are consistent
        sz_check = (sz_dims == sz_obj);
        sz_check(dim) = []; % obviously for catting dim, sizes won't be consistent
        if(~all(sz_check))
            error('rome_var:cat:InconsistentDimensions', 'Object sizes must be consistent');
        end            
    end
    
    % compute number of vars (mapped and unmapped)
    num_mapped_var_arr(ii) = cur_obj.NumMappedVars;
    num_unmapped_var_arr(ii) = cur_obj.NumUnmappedVars;
    
    % do the same for uncertain variables
    num_mapped_rand_var_arr(ii) = cur_obj.NumMappedRandVars;
    num_unmapped_rand_var_arr(ii) = cur_obj.NumUnmappedRandVars;
    
    % flags
    is_pure_rand_arr(ii) = cur_obj.IsRand;
    is_pure_certain_arr(ii) = cur_obj.IsCertain;
end

% preprocess
num_unmapped_var_arr(logical(is_pure_rand_arr)) = max(num_unmapped_var_arr);
num_unmapped_rand_var_arr(logical(is_pure_certain_arr)) = max(num_unmapped_rand_var_arr);

% allocate memory for new variable
sz_dims(dim) = sum(cat_dim_arr);
out_var_obj = rome_var(sz_dims);
out_var_obj.NumUnmappedVars = min(num_unmapped_var_arr);
out_var_obj.NumUnmappedRandVars = min(num_unmapped_rand_var_arr);

% compute total number of vars
num_total_var_arr       = num_mapped_var_arr + num_unmapped_var_arr;
num_total_rand_var_arr  = num_mapped_rand_var_arr + num_unmapped_rand_var_arr;

% number of total variables in output (mapped + unmapped)
max_total_vars      = max(num_total_var_arr);
max_total_rand_vars = max(num_total_rand_var_arr);

% number of mapped variables in output object
num_mapped_vars               = max_total_vars - out_var_obj.NumUnmappedVars;
out_var_obj.NumMappedRandVars = max_total_rand_vars - out_var_obj.NumUnmappedRandVars;

% define pre and post zero padding array
num_pre_zero_pad = num_unmapped_var_arr - out_var_obj.NumUnmappedVars;
num_post_zero_pad= max_total_vars - num_total_var_arr;
num_pre_zero_pad_rand = num_unmapped_rand_var_arr - out_var_obj.NumUnmappedRandVars;
% num_post_zero_pad_rand= max_total_rand_vars - num_total_rand_var_arr;

% construct permutation indices
inner_sz = prod(sz_dims(1:dim-1));
outer_sz = prod(sz_dims(dim+1:end));
perm_ind = zeros(inner_sz * sz_dims(dim), outer_sz);

row_start_ind = 1 + inner_sz * cumsum([0, cat_dim_arr(1:end-1)]);
row_end_ind   = inner_sz * cumsum(cat_dim_arr);

actual_start_ind  = 1 + (inner_sz * outer_sz) * [0, cumsum(cat_dim_arr(1:end-1))];
actual_end_ind    = (inner_sz * outer_sz) * cumsum(cat_dim_arr);

% construct the final output
out_var_obj.BiAffineMap = [];
for ii = 1:num_objs
    % ignore empty matrices
    cur_obj = varargin{ii};
    if(isempty(cur_obj))
        continue;
    end
    
    cur_num_out_total_vars = cur_obj.TotalSize;
    
    % allocate memory
    insert_map = spalloc(cur_num_out_total_vars, (num_mapped_vars + 1) * (out_var_obj.NumMappedRandVars + 1) , nnz(cur_obj.BiAffineMap));
    
    % inner indices
    inner_sz  = 1 + cur_obj.NumMappedVars + num_pre_zero_pad(ii) + num_post_zero_pad(ii);
    inner_ind = [1, num_pre_zero_pad(ii) + 1 + (1:cur_obj.NumMappedVars)];     
    
    % make indices
    inner_ind_expand = outer_sum(inner_ind', 0:inner_sz:(inner_sz*(cur_obj.NumMappedRandVars-1)));
        
    % make the complete index vector
    ind_vec = [inner_ind, inner_sz * (num_pre_zero_pad_rand(ii) + 1) + inner_ind_expand(:)'];
                 
    % assign
    insert_map(:, ind_vec) = cur_obj.BiAffineMap;    
    out_var_obj.BiAffineMap = [out_var_obj.BiAffineMap; insert_map];
    
    if(~isempty(perm_ind))
        perm_ind(row_start_ind(ii):row_end_ind(ii), :) = ...
            reshape(actual_start_ind(ii):actual_end_ind(ii), ...
            row_end_ind(ii) - row_start_ind(ii) + 1, outer_sz);
            
    end
end

% do permutation
out_var_obj.BiAffineMap = out_var_obj.BiAffineMap(perm_ind(:), :);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
