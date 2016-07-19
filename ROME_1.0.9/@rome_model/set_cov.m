% function set_cov(model_obj, nprior_vars, cov_mat, cov_mix)
function set_cov(model_obj, var_obj, cov_mat)

% ROME_MODEL\SET_COV Internal method that supplies covariance information
% about a uncertain variable
%
% This runs for global_x(nprior_vars+1:nprior_vars+length(mean_val)).
% notice that cov_val is to be a symmetric matrix, preferably sparse. The
% calling function should convert all other types of input (scalar /
% vector) into matrix form before calling this function.
%
% Modification History: 
% 1. Joel 

% % compute number of previous covar results
% num_prev_cov = size(model_obj.rndCovar, 1);
% num_unknown_cov  = pos(nprior_vars - num_prev_cov);
% 
% % MATLAB auto-allocates
% ind = {nprior_vars+(1:size(cov_mat, 1)), nprior_vars+(1:size(cov_mat, 2))};
% model_obj.rndCovar(ind{:}) = cov_mat;
% 
% model_obj.rndCovar(num_prev_cov + (1:num_unknown_cov), num_prev_cov + (1:num_unknown_cov)) ...
%     = spdiags(NaN(num_unknown_cov, 1), 0, num_unknown_cov, num_unknown_cov);
 

model_obj.rndCovar = cov_mat;
model_obj.rndCovarMix = var_obj;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
