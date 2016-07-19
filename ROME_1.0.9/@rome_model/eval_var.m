function val = eval_var(model_obj, var_obj)

% ROME_MODEL\EVAL_VAR Evaluates the current value of a variable using the
% solution found from the model
%
% Modified by:
% Joel

if(isa(var_obj, 'rome_var'))
    % Find out range of the var_obj
    val_index = var_obj.NumUnmappedVars + (1:var_obj.NumMappedVars);
    
    if(~isempty(val_index))
        % CASE: Pure Certain or LDR Vars
        % extract relevant xsol values and substitute into variable
        val = var_obj.insert(model_obj.XSol(val_index));
    else
        % CASE: Pure Uncertain Vars
        % TEST
        val = var_obj.BiAffineMap;        
    end
    
    % return as a full matrix
    val = full(val);
    
elseif(isa(var_obj, 'rome_constraint_index'))
    
    % find non-inf bounds
    num_primal_LC = model_obj.LC.TotalSize;
    num_LB = length(~isinf(model_obj.LB));
    
    % modify index based on constraint type
    for ii = 1:numel(var_obj)
        if(var_obj(ii).Type == rome_constraint_index.LB)
            var_obj(ii).Index = var_obj(ii).Index + num_primal_LC;
        elseif((var_obj(ii).Type == rome_constraint_index.UB))
            var_obj(ii).Index = var_obj(ii).Index + num_primal_LC + num_LB;
        end
    end
        
    % vertical concatenation of all indices
    ind_arr = vertcat(var_obj.Index);
    
%     % modify index based on type
%     type_arr = vertcat(var_obj.Type);
%     LB_ind = (type_arr == rome_constraint_index.LB);
%     ind_arr(LB_ind) = ind_arr(LB_ind) + num_primal_LC;
%     
%     UB_ind = (type_arr == rome_constraint_index.UB);    
%     ind_arr(UB_ind) = ind_arr(UB_ind) + num_primal_LC + num_LB;    
    
    % extract from dual variable
    val = model_obj.DualVars(ind_arr);
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
