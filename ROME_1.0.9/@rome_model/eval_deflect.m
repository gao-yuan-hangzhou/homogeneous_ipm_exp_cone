function [nonlin_coeffs, nonlin_values] = eval_deflect(model_obj, var_obj)

% ROME_MODEL\EVAL_DEFLECT Evaluates the deflection of an object from the
% solution found in the model
%
% USAGE:
% [nonlin_coeffs, nonlin_values] = eval_deflect(model_obj, var_obj);
%
%
% Modified by:
% 1. Joel (14 Apr 2009)
%

% initialize outputs
nonlin_values = [];
nonlin_coeffs = [];

% Quick check if there were any successful deflections. 
% If so, compute their value, otherwise, return immediately

if(isempty(model_obj.ldrRemLC))
    return;
else
    obj = model_obj.ldrRemLC;  
    
    % notice transpose
    nonlin_values = squeeze(obj.insert(model_obj.XSol(obj.NumUnmappedVars + 1:obj.NumMappedVars)))';
    nonlin_values = nonlin_values(:)';
    
    % reshape nonlin values to have N+1 cols
    N = obj.NumMappedRandVars;
    M = length(nonlin_values) ./ (N+1);
    nonlin_values = reshape(nonlin_values', N + 1, M)';
    
    % Fix for uncertain means (27 Mar 2010)
    ind_not_mean = [true; ~model_obj.ZIsMean];
    nonlin_values = nonlin_values(:, ind_not_mean);
end
    
if(~isa(var_obj, 'rome_var'))    
    error('rome_model:eval_deflect:OnlyAcceptRomeVar', 'var_obj must be of type rome_var');
else
    if(~var_obj.IsLDR)
        return;
        %error('rome_model:eval_deflect:OnlyDeflectLDR', 'Only LDR variables have deflections');
    end
    
    % get the certain part of the var_object only
    certain_var_obj = strip_rand(var_obj);
    
    % extract relevant cols from the input object
    var_cols = model_obj.ldrDefInd - certain_var_obj.NumUnmappedVars;
    capped_cols = var_cols(var_cols >= 1 & var_cols <= certain_var_obj.NumMappedVars + 1);
    
    var_coeffs = [zeros(certain_var_obj.TotalSize, sum(var_cols < 1)), ...
                  full(certain_var_obj.BiAffineMap(:, capped_cols)), ...
                  zeros(certain_var_obj.TotalSize, sum(var_cols > certain_var_obj.NumMappedVars + 1))];
    
    % notice use of matrix multiply.
    nonlin_coeffs = var_coeffs * model_obj.ldrDefCoeff;
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
