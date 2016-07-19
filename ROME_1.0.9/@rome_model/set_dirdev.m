function set_dirdev(model_obj, var_obj, dirdev, dir_flag)

% ROME_MODEL\SET_DIRDEV Internal method that supplies directional deviation
% information about a uncertain variable. 
%
% If dir_flag > 0, sets the forward deviation, 
% If dir_flag < 0, sets the backward deviation, 
% If dir_flag = 0, yields an error
%
% This runs for global_x(nprior_vars+1:nprior_vars+length(dirdev)).
% notice that dirdev is to be a vector, and the calling function should
% convert all scalar mean_vals into vectors
%
% Modification History: 
% 1. Joel 

if(dir_flag > 0)
    model_obj.rndFDev = dirdev;
    model_obj.rndDirDevMix = var_obj;
elseif(dir_flag < 0)
    model_obj.rndBDev = dirdev;
    model_obj.rndDirDevMix = var_obj;
else
    error('rome_model:set_dirdev:UnknownDirFlag', 'Direction Flag cannot be zero.');
end

% OLD CODE
% if(dir_flag > 0)
%     % compute number of forward deviations
%     num_prev_pdev = length(model_obj.rndFDev);
%     num_unknown_pdev  = pos(nprior_vars - num_prev_pdev);
% 
%     % MATLAB auto-allocates
%     model_obj.rndFDev(nprior_vars + (1:length(dirdev)), 1) = dirdev;
%     model_obj.rndFDev(num_prev_pdev + (1:num_unknown_pdev)) = NaN;
% elseif(dir_flag < 0)
%     % compute number of backward deviations
%     num_prev_qdev = length(model_obj.rndBDev);
%     num_unknown_qdev  = pos(nprior_vars - num_prev_qdev);
% 
%     % MATLAB auto-allocates
%     model_obj.rndBDev(nprior_vars + (1:length(dirdev)), 1) = dirdev;
%     model_obj.rndBDev(num_prev_qdev + (1:num_unknown_qdev)) = NaN;
% else
%     error('rome_model:set_dirdev:UnknownDirFlag', 'Direction Flag cannot be zero.');
% end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
