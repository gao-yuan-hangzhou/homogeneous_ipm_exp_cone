function add_new_var(model_obj, var_obj)

% PROFMODEL\ADD_NEW_VAR Appends a new variable to the current model
%
%
% Modification History:
% 1. Joel

if(var_obj.IsRand())
    % append rndZ with zeros
    model_obj.ZRnd = [model_obj.ZRnd; zeros(var_obj.TotalSize, 1)];   % 2 times because the next one is for mean

    % append Vartype with Character flag
    switch(var_obj.Continuity)
        case rome_constants.CONTINUOUS
            cur_vartype = 'C';
        otherwise
            error('rome_model:add_new_var:UnrecognizedVarType', ...
                'Uncertain Variable type must be CONTINUOUS');
    end

    % assign variable type
    model_obj.rndVarType = [model_obj.rndVarType; repmat(cur_vartype, var_obj.TotalSize, 1)];
  
else
    % append x_sol with zeros
    model_obj.XSol = [model_obj.XSol; zeros(var_obj.TotalSize, 1)];

    % append Vartype with Character flag
    switch(var_obj.Continuity)
        case rome_constants.CONTINUOUS
            cur_vartype = 'C';
        case rome_constants.INTEGER
            cur_vartype = 'I';
        case rome_constants.BINARY
            cur_vartype = 'B';
        otherwise
            error('rome_model:add_new_var:UnrecognizedVarType', ...
                'Variable type must be one of CONTINUOUS, INTEGER, or BINARY');
    end

    % assign variable type
    model_obj.VarType = [model_obj.VarType; repmat(cur_vartype, var_obj.TotalSize, 1)];
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
