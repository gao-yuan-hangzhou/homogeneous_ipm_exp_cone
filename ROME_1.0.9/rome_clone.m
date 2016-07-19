function out_obj = rome_clone(obj)

% ROME_CLONE Makes a clone of the object having the same dependency
% structure.
%
%
% Input  : obj (should be an LDR rome_var)
% Returns: out_obj, a newly-created variable with the same dependency
% structure as obj.
%
% Modification History: 
% 1. Joel (11 Feb 2009)

% check for pathological case
if(obj.IsCertain)
    out_obj = rome_model_var(size(obj));
    return;
elseif(obj.IsRand)
    out_obj = rome_rand_model_var(size(obj));
    return;
end

% get the current model
h = rome_get_current_model();

% get the primitive uncertain variables for obj
z = h.get_rand_vars(obj);

% find dependency structure by building info-index table
inner_step = obj.NumMappedVars + 1;
[r, c] = find(obj.BiAffineMap);

% notice vectorize of r and c
unique_ind = unique([r(:), ceil(c(:) ./ inner_step)], 'rows');
info_index_table = logical(accumarray(unique_ind, 1));

% expand the information index table
info_index_table = [info_index_table, false(obj.TotalSize, pos(z.TotalSize - size(info_index_table,2) + 1))];

% (June 26 2012) Edited out in ROME 1.0.9: Why did we do this in the first place? Joel to check
% expand_info_index_table = repmat(info_index_table(:), 1, out_obj.NumMappedVars + 1);
% expand_info_index_table = blk_transpose(expand_info_index_table, obj.TotalSize, out_obj.NumMappedVars + 1);

% % zero out non-dependent variables
% out_obj.BiAffineMap(~expand_info_index_table) = 0;

% make a copy of the object
out_obj = rome_linearrule(size(obj), z, 'Pattern', info_index_table);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
