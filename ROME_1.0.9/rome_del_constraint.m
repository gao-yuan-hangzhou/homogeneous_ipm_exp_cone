function rome_del_constraint(constraint_index_obj)

% ROME_DEL_CONSTRAINT Deletes an existing constraint in the moodel
%
% Usage:
%   mod_constraint(ROME_CONSTRAINT_INDEX_OBJ)
% 
% Example:
%        m = rome_model('Test Model);
%       ...
%        p = rome_constraint(x + 3 * y <= 15);
%       ...
%       
%        rome_del_constraint(p);    <-- removes Constraint p from the model
%
% Notes:
% Modification History: 
% 1. Joel 

% error checking
if(~isa(constraint_index_obj, 'rome_constraint_index'))
    error('rome_model:mod_constraint:InvalidConstraintIndex', 'Number of Prior Vars must be a non-negative integer.');
end

% get the current model 
h_curr_model = rome_get_current_model();

% switch based on constraint type
if(constraint_index_obj.Type == rome_constraint_index.LC)
    % find indices which are NOT in set to be deleted
    remain_ind = (1:h_curr_model.LC.BiAffineMap.TotalSize).';
    remain_ind = setdiff(remain_ind, constraint_index_obj.Index);
    
    % change the affineMap for the Linear constraint
    h_curr_model.LC.BiAffineMap = h_curr_model.LC.BiAffineMap(remain_ind, :);
    
elseif(constraint_index_obj.Type == rome_constraint_index.LB)
    % By this juncture, the LB array should have already been allocated
    % Removing a LB constraint is same as making x >= -inf
    h_curr_model.LB(constraint_index_obj.Index) = -Inf;    
    
elseif(constraint_index_obj.Type == rome_constraint_index.UB)
    % By this juncture, the LB array should have already been allocated
    % Removing a LB constraint is same as making x <= inf
    h_curr_model.UB(constraint_index_obj.Index) = Inf;
else
    error('TEMP', 'Constraint Index Type not recognized');
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
