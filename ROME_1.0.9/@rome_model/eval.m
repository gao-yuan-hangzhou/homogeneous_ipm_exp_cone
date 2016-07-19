function sol_obj = eval(model_obj, var_obj)

% ROME_MODEL\EVAL Evaluates the current full variable with deflected
% components. Returns a rome_sol object
%
% Modified by:
% Joel

if(isa(var_obj, 'rome_var'))
    % catch simple case when no deflection
    if(var_obj.IsCertain)
        sol_obj = model_obj.eval_var(var_obj);
        return;
    end    
    
    % standard case with uncertainty
    sol_obj = rome_sol;
    [sol_obj.DeflectCoeff, sol_obj.DeflectAffineMap] = model_obj.eval_deflect(var_obj);
    sol_obj.LDRAffineMap = squeeze(model_obj.eval_var(var_obj(:)));
    sol_obj.Size = var_obj.Size;
    
    % catch abberant case
    if(all(var_obj.TotalSize == 1))
        if(size(sol_obj.LDRAffineMap, 1) ~= 1)
            sol_obj.LDRAffineMap = sol_obj.LDRAffineMap';
        end
    end
    
    % Mar 19 2010: add trailing zeros to handle non-full uncertainty
    % dependence
    sol_obj.LDRAffineMap = [sol_obj.LDRAffineMap, ...
        zeros(size(sol_obj.LDRAffineMap, 1), model_obj.NumRandVars - sum(model_obj.ZIsMean) - size(sol_obj.LDRAffineMap, 2) + 1)];
   
else
    error('rome_model:eval_all:OnlyRomeVar', 'Only accepts rome_var objects.');
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
