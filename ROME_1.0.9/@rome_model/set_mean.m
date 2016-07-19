function set_mean(model_obj, nprior_vars, mean_val)

% ROME_MODEL\SET_MEAN Internal method that sets mean value information
% about a uncertain variable
%
% This runs for global_x(nprior_vars+1:nprior_vars+length(mean_val)).
% notice that mean_val is to be a vector, and the calling function should
% convert all scalar mean_vals into vectors
%
% Modification History: 
% 1. Joel 

% compute number of previous means and unknown means
num_prev_means = length(model_obj.rndMean);
num_unknown_means  = pos(nprior_vars - num_prev_means);

% % notice that we use use (len, 1) to ensure mean is a col vector
% if(isnumeric(mean_val))
%     % MATLAB auto-allocates
%     model_obj.rndMean(nprior_vars + (1:length(mean_val)), 1) = mean_val; 
%     model_obj.rndUnknownMeanInd(nprior_vars + (1:length(mean_val)), 1) = false;
% else
%     model_obj.rndMean(nprior_vars + (1:length(mean_val)), 1) = ...
%         mean_val.NumUnmappedRandVars + (1:mean_val.NumMappedRandVars);
%     model_obj.rndUnknownMeanInd(nprior_vars + (1:length(mean_val)), 1) = true;
% end
% model_obj.rndMean(num_prev_means + (1:num_unknown_means), 1) = NaN;


% TEMP!
% notice that we use use (len, 1) to ensure mean is a col vector
model_obj.rndMean(nprior_vars + (1:length(mean_val)), 1) = mean_val;

if(num_unknown_means > 0)
    error('rome_model:set_mean:OutOfBounds', 'Cannot set means out of bounds');
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
