function mean_val = get_mean(model_obj, var_obj)

% PROF_MODEL\GET_MEAN Internal method that retrieves means of uncertain
% variables 
%
% This runs for global_x(nprior_vars+1:nprior_vars+len). Returns NaN if out-of bounds.
%
% Modification History: 
% 1. Joel 

% find start and end indices
start_index = var_obj.NumUnmappedRandVars + 1;
end_index   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;

% % compute number of NaNs
% num_NaNs = pos(end_index - numel(model_obj.rndMean));
% end_index = end_index - num_NaNs;   % need to shorten end_index accordingly if number of NaNs > 0

% % check if all mean values are deterministic
% ind = logical(model_obj.rndUnknownMeanInd(start_index:end_index));
% if(any(ind))
%     mean_val = rome_empty_var(end_index - start_index + 1 + num_NaNs);
%     mean_val(find(ind)) = model_obj.get_rand_vars(model_obj.rndMean(ind));
%     mean_val(find(~ind)) = model_obj.rndMean(~ind);    
% else
%     % All Deterministic - just return
%     mean_val = [model_obj.rndMean(start_index:end_index); NaN(num_NaNs, 1)];
% end

% TEMP!
if(pos(end_index - model_obj.rndMean.TotalSize) > 0)
    error('rome_model:get_mean:OutOfBounds', 'Trying to get mean that is out of bounds');
end

% grab mean vals
mean_val = model_obj.rndMean(start_index:end_index);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
