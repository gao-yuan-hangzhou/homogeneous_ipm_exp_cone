 function obj = repmat(obj, pattern)

% ROME_VAR\REPMAT Repeats a rome_var object in a specified pattern
%
% new_obj = reshape(obj, sz)
% new_obj = reshape(obj, N1, N2, N3... )
%
% Modification History: 
% 1. Joel 

% % parse arguments
% if(numel(varargin) == 1)
%     new_sz = varargin{1};
% else
%     new_sz = horzcat(varargin{:});
% end

% error check
if(any(int64(pattern) ~= pattern))
    error('rome_var:reshape:ReqIntegerSize', 'Size must be real integer values');
end

% vectorize pattern
pattern  = pattern(:)';

% squeeze out trailing singletons in pattern
while((length(pattern) > 2) && (pattern(end) == 1))
    pattern = pattern(1:end-1);
end

% % METHOD 2
% % expand size to be same dimensions as pattern or vice versa
% old_sz = size(obj);
% new_sz = old_sz;
% diff_dim = numel(pattern) - numel(old_sz);
% 
% % equalize dimensions of size
% if(diff_dim > 0)
%     new_sz = [new_sz, ones(1, diff_dim)];
% else
%     pattern = [pattern, ones(1, -diff_dim)];
% end
% new_sz = new_sz .* pattern;

% % construct permutation indices
% perm_ind = repmat(reshape(1:obj.TotalSize, old_sz), pattern);
%  
% % do repeating
% new_map = obj.BiAffineMap;
%  
% % permute
% obj.BiAffineMap = new_map(perm_ind(:), :);
%  
% % if all ok, change the Size of the object
% obj.Size = new_sz;

% METHOD 0
% expand size to be same dimensions as pattern or vice versa
old_sz = size(obj);
diff_dim = numel(pattern) - numel(old_sz);

% equalize dimensions of size
if(diff_dim > 0)
    old_sz = [old_sz, ones(1, diff_dim)];
else
    pattern = [pattern, ones(1, -diff_dim)];
end

% do repeating OLD (super fast)
ndims    = numel(pattern);
old_cols = size(obj.BiAffineMap, 2);

% define vector of rows to grab for repeating
grab_row_vec = [old_sz; pattern];
grab_row_vec = reshape(cumprod(grab_row_vec(:)), 2, ndims);

% iterate through each dimension and expand the map accordingly
new_map  = obj.BiAffineMap;
for ii = 1:ndims
    new_map = blk_transpose(new_map, grab_row_vec(1, ii), old_cols);
    new_map = repmat(new_map, [pattern(ii), 1]);
    new_map = blk_transpose(new_map, grab_row_vec(2, ii), old_cols);
end

% assign map and size
obj.BiAffineMap = new_map;
obj.Size = old_sz .* pattern;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
