function rome_minimize(obj)

% ROME_MINIMIZE Sets the objective in the current model
%
%   obj must be a scalar-valued rome_var or a real scalar.
%
%   This should be an internal function called by the engine. 
%
% Modification History: 
% 1. Joel 

% error_check
if(~isscalar(obj))
    error('rome_minimize:InvalidArg', 'Input must be scalar-valued.');
end

% Get handle to the current model
h_curr_model = rome_get_current_model();

% check that object is not LDR. If is, convert to epigraph 
if(obj.IsLDR)
    t = rome_model_var;
    rome_constraint(obj <= t);
    h_curr_model.ObjFn = t;
else
    h_curr_model.ObjFn = obj;
end

h_curr_model.MinMaxFlag = rome_model.MINIMIZE;



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
