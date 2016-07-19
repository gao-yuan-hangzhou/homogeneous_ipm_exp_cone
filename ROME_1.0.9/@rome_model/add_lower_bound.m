function add_lower_bound(obj, nprior_vars, b, is_rand_var, nskip)

% ROME_MODEL\ADD_LOWER_BOUND Internal method that adds a lower bound
% constraint of the form x >= b.
%
% If is_rand_var is 1, will treat the bound as a uncertain variable,
% otherwise, if is_rand_var is 0 or omitted, will treat as regular variable
%
% The constraint runs for global_x(nprior_vars+1:nskip:nprior_vars+length(b)).
%
% This function should also handle redundant constraints. 
% 
% E.g. 
%       rome_constraint(x >= 4) 
%       rome_constraint(x >= 5) <-- will only retain this one 
%
% Modification History: 
% 1. Joel 

if(nargin < 5)
    nskip = 1;          % by default, skip = 1
end
if(nargin < 4)
    is_rand_var = 0;    % by default, not a uncertain variable
end

if(~isscalar(nprior_vars) || ~isnumeric(nprior_vars) || nprior_vars < 0)
    error('rome_model:add_lower_bound:InvalidPriorVars', 'Number of Prior Vars must be a non-negative integer.');
end
if(~isvector(b))
    error('rome_model:add_lower_bound:InvalidBound', 'Bounding vector must be a vector');
end

end_index = nprior_vars + nskip * length(b);

if(~is_rand_var)
    % Case 1: When its a normal variable
    num_append = end_index - length(obj.LB);

    % expand array if necessary
    if(num_append > 0)
        obj.LB = [obj.LB; repmat(-Inf, num_append, 1)];
    end

    % handle redundant constraints
    cur_LB_sel = obj.LB((nprior_vars+1):nskip:end_index);
    if(~isinf(cur_LB_sel))
        warning('rome_model:add_lower_bound:RedundantConstraint', ...
            'Redundant diagonal constraint: dual variable will be with respect to tighter constraint.');
    end
    new_active_indices = find(b > cur_LB_sel);
    obj.LB(nprior_vars + nskip * (new_active_indices-1) + 1) = b(new_active_indices);
else    
    % Case 2: When its a uncertain variable
    num_append = end_index - length(obj.rndLB);

    % expand array if necessary
    if(num_append > 0)
        obj.rndLB = [obj.rndLB; repmat(-Inf, num_append, 1)];
    end

    % handle redundant constraints
    cur_LB_sel = obj.rndLB((nprior_vars+1):end_index);
    if(~isinf(cur_LB_sel))
        warning('rome_model:add_lower_bound:RedundantConstraint', ...
            'Redundant diagonal constraint: dual variable will be with respect to tighter constraint.');
    end
    new_active_indices = find(b > cur_LB_sel);
    obj.rndLB(nprior_vars + new_active_indices) = b(new_active_indices); 
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
