function rome_mod_constraint(constraint_index_obj, new_constraint)

% ROME_MOD_CONSTRAINT Modifies an existing constraint for the CURRENT
% model. 
%
% Usage:
%   mod_constraint(ROME_CONSTRAINT_INDEX_OBJ, NEW_ROME_VAR_CONSTRAINT)
% 
% Example:
%        m = rome_model('Test Model);
%       ...
%        p = rome_constraint(x + 3 * y <= 15);
%       ...
%       
%        mod_constraint(p, x + 4 * y <= 4);
%
% Notes:
%   NEW_ROME_VAR_CONSTRAINT must be a valid rome_var object, and also should have a
%   proper conic restriction. A simple way to do this is to use overloaded
%   <=, >= or even == to constrain a rome_var object.
%
% Modification History: 
% 1. Joel 

% error checking
if(~isa(constraint_index_obj, 'rome_constraint_index'))
    error('rome_model:mod_constraint:InvalidConstraintIndex', 'Number of Prior Vars must be a non-negative integer.');
end
if(~isa(new_constraint, 'rome_var'))
    error('rome_model:mod_constraint:InvalidNewConstraint', 'New constraint must be a rome_var object');
end
if(new_constraint.Cone == rome_constants.NO_CONE)
    warning('rome_model:mod_constraint:NoEffectiveConstraint', 'New constraint is a free rome_constants. Is this what you wanted?');
end

% get the current model 
h_curr_model = rome_get_current_model();

% assume new constraint is a Linear Constraint
if(constraint_index_obj.Type ~= rome_constraint_index.LC)
    error('TEMP');
else
    % notice that this will throw an error if the sizes are not matched
    h_curr_model.LC.BiAffineMap(constraint_index_obj.Index, :) = new_constraint.BiAffineMap;
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
