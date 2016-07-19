function [cov_mat, cov_mix] = get_cov(model_obj, var_obj)

% ROME_MODEL\GET_COV Internal method that retrieves covariance of uncertain
% variables 
%
%   var_obj should be a uncertain rome_var, indicating the desired range. 
%   var_obj can also be a 2x1 vector = [start_index, end_index]
%
% Modification History: 
% 1. Joel 

% % find start and end indices
% start_index = var_obj.NumUnmappedRandVars + 1;
% end_index   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;
% 
% % compute number of NaNs
% num_NaNs = pos(end_index - size(model_obj.rndCovar, 1));
% end_index = end_index - num_NaNs;   % need to shorten end_index accordingly if number of NaNs > 0
% 
% % return covariance
% cov_mat = model_obj.rndCovar(start_index:end_index, start_index:end_index);
% cov_mat(end_index+1:num_NaNs, end_index+1:num_NaNs) = spdiags(NaN(num_NaNs, 1), 0, num_NaNs, num_NaNs);
% 
% % return covariance mixing matrix
% cov_mix = model_obj.rndCovarMix(start_index:end_index, start_index:end_index);

mix_obj = model_obj.rndCovarMix;
if(nargin > 1)
    if(isa(var_obj, 'rome_var'))
        start_index = var_obj.NumUnmappedRandVars + 1;
        end_index   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;
    else
        start_index = var_obj(1);
        end_index = var_obj(2);
    end
    
    % check
    if((mix_obj.NumUnmappedRandVars < start_index - 1) ||  (mix_obj.NumUnmappedRandVars + mix_obj.NumMappedRandVars > end_index))
        error('rome_model:get_cov', 'Indices of covariance to be extracted must not be a subset of supplied covariance!');
    end
    
    % check
    if(end_index > model_obj.NumRandVars)
        error('rome_model:get_cov', 'end_index is larger than number of uncertain variables');
    end
else
    start_index = mix_obj.NumUnmappedRandVars + 1;
    end_index   = mix_obj.NumUnmappedRandVars + mix_obj.NumMappedRandVars;
end

% get the covariance and mixing matrix from the model object
cov_mat = model_obj.rndCovar;
cov_mix = mix_obj.BiAffineMap(:, 2:end);

% zero-pad the cov_mixing matrix
n = size(cov_mix, 1);
m = mix_obj.NumUnmappedRandVars - start_index + 1;
cov_mix = [zeros(n, m), cov_mix, zeros(n, end_index - mix_obj.NumUnmappedRandVars - mix_obj.NumMappedRandVars)];

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
