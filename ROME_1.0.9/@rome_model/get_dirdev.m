function [dirdev_val, dirdev_mix] = get_dirdev(model_obj, var_obj, dir_flag)

% ROME_MODEL\GET_DIRDEV Internal method that retrieves directional deviations 
% of uncertain variables 
% 
% If dir_flag > 0, sets the forward deviation, 
% If dir_flag < 0, sets the backward deviation, 
% If dir_flag = 0, yields an error
%
%
% Modification History: 
% 1. Joel 

mix_obj = model_obj.rndDirDevMix;
if(nargin > 1)
    if(isa(var_obj, 'rome_var'))
        start_index = var_obj.NumUnmappedRandVars + 1;
        end_index   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;
    else
        start_index = var_obj(1);
        end_index = var_obj(2);
    end
    
    % check need
    if((mix_obj.NumUnmappedRandVars < start_index - 1) ||  (mix_obj.NumUnmappedRandVars + mix_obj.NumMappedRandVars > end_index))
        error('rome_model:get_dirdev:MustGetAllCovar', ...
                'Indices of dir deviation to be extracted must not be a subset of previously input dir deviation!');
    end
    
    % check
    if(end_index > model_obj.NumRandVars)
        error('rome_model:get_dirdev:OutOfBounds', 'end_index is larger than number of uncertain variables');
    end
else
    start_index = mix_obj.NumUnmappedRandVars + 1;
    end_index   = mix_obj.NumUnmappedRandVars + mix_obj.NumMappedRandVars;
end

% get the covariance and mixing matrix from the model object
if(dir_flag > 0)
    dirdev_val = model_obj.rndFDev;
elseif(dir_flag < 0)
    dirdev_val = model_obj.rndBDev;
else
    error('rome_model:get_dirdev:UnknownDirFlag', 'Direction Flag cannot be zero.');
end
dirdev_mix = mix_obj.BiAffineMap(:, 2:end);   

% zero-pad the mixing matrix
n = size(dirdev_mix, 1);
m = mix_obj.NumUnmappedRandVars - start_index + 1;
dirdev_mix = [zeros(n, m), dirdev_mix, zeros(n, end_index - mix_obj.NumUnmappedRandVars - mix_obj.NumMappedRandVars)];

% put in the linear component
dirdev_mix = [mix_obj.BiAffineMap(:, 1), dirdev_mix];

% OLD CODE
% % find start and end indices
% start_index = var_obj.NumUnmappedRandVars + 1;
% end_index   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;
% 
% if(dir_flag > 0)
%     % compute number of NaNs
%     num_NaNs = pos(end_index - numel(model_obj.rndFDev));
%     end_index = end_index - num_NaNs;   % need to shorten end_index accordingly if number of NaNs > 0
%     
%     % return directional deviation value
%     dirdev_val = [model_obj.rndFDev(start_index:end_index); NaN(num_NaNs, 1)];
%     
% elseif(dir_flag < 0)
%     % compute number of NaNs
%     num_NaNs = pos(end_index - numel(model_obj.rndBDev));    
%     end_index = end_index - num_NaNs;   % need to shorten end_index accordingly if number of NaNs > 0    
%     
%     % return directional deviation value
%     dirdev_val = [model_obj.rndBDev(start_index:end_index); NaN(num_NaNs, 1)];
% else
%     error('rome_model:get_dirdev:UnknownDirFlag', 'Direction Flag cannot be zero.');
% end





% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
