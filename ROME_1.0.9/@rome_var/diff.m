function out_var_obj = diff(A, N, dim)

% ROME_VAR\DIFF Difference of rome_vars
%
%   C = diff(A, N, dim)
%   C = diff(A, N);
%   C = diff(A)
%
% Modification History: 
% 1. Joel 

% catch pathological case
if(isempty(A))
    out_var_obj = [];
    return;
end

% default arguments
if(nargin < 3) 
    dim = find(size(A) > 1);
    if(~isempty(dim))
        dim = dim(1);
    else
        dim = 1;
    end
end
if(nargin < 2)
    N = [];
end

if((~isempty(N)) && (N ~= 1))
    error('rome_var:diff:NotImplemented', 'Haven''t implemented higher order differencing yet');
end

% get original size
orig_size = size(A);

if(dim < 1 || ~isscalar(dim) || floor(dim) ~= dim)
    error('rome_var:diff:InvalidArg', ...
        'Differencing dimension must be positive scalar integer within indexing range');
elseif(dim > length(orig_size) || orig_size(dim) == 1)
    % differencing over a larger dimension or a singleton dimension returns null matrix
    out_var_obj = [];
    return;
end

% create new size vector, squeeze out trailing singleton dimensions
new_size = orig_size;
new_size(dim) = orig_size(dim) - 1;     % the size of the differencing dimension reduces by 1
if(new_size(end) == 1)
    new_size = new_size(1:end-1);
end

% allocate memory for new object
out_var_obj = rome_var(new_size);
out_var_obj.NumUnmappedVars = A.NumUnmappedVars; % no change to unmapped variables
out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars; % no change to unmapped variables
out_var_obj.NumMappedRandVars = A.NumMappedRandVars; % no change to unmapped variables

% create premultiplying matrix
inner_sz = prod(orig_size(1:dim-1));
tot_vars = prod(orig_size);
premult_mat = spdiags([-ones(tot_vars, 1), ones(tot_vars, 1)], [0, inner_sz], tot_vars, tot_vars);

% remove some rows from premult_mat
del_ind = outer_sum(((inner_sz * (orig_size(dim) - 1)):inner_sz*orig_size(dim):tot_vars)', 1:inner_sz);
premult_mat(del_ind(:), :) = [];

% compute
out_var_obj.BiAffineMap = premult_mat * A.BiAffineMap;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
