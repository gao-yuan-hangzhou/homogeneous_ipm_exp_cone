function [bound_obj] = rome_pwl_covar_bound(y, b_vec, a_vec)

% ROME_PWL_COVAR_BOUND Generates the mean-covariance bound for a piece-wise linear 
% function of an LDR. (TESTING)
%
% sup_{z_mean}{0.5 (y_0 + y'z_mean) 
%             + 0.5 * sqrt( (y_0 + y'z_mean)^2 + y' * Covar * y ) }
%
% coeff should be formatted as 
%   coeff = [b_1, a_1; b_2, a_2; b_3, a_3 ... ; b_k, a_k] 
%
% Input  : y (should be an LDR rome_var)
% Returns: Bounding variable
%
% Modification History: 
% 1. Joel (13 Feb 2009)

% get the current model
h = rome_get_current_model();

% make variables
w = rome_model_var(y.TotalSize);
s = rome_model_var(y.TotalSize);
t = rome_model_var(y.TotalSize);
z = rome_model_var(y.TotalSize, 'Cone', rome_constants.NNOC);
bound_obj = w + s;

% apply the k constraints
for ii = 1:size(b_vec, 2)
    if(ndims(a_vec) > 2)
        rome_constraint(w - a_vec(:, :, ii).^2 * z - a_vec(:, :, ii) * (t + y.mean) - b_vec(:, ii) >= 0);
    else
        % optimization to handle case of elements of a_vec being scalar
        rome_constraint(w - a_vec(ii).^2 * z - a_vec(ii) * (t + y.mean) - b_vec(:, ii) >= 0);
    end
end

% get the factored covariance matrix
[proj_covar, F_mix] = h.get_cov(y);
z_std = chol(proj_covar); % projected covariance
% z_std = full(proj_covar)^(0.5); % projected covariance % HACK!

% apply the SOC constraint
X = rome_model_var(size(F_mix, 1), y.TotalSize);
rome_constraint(strip_certain(y) == F_mix' * X); 

tmp = [(s + z)'; z_std * X; t'; (s - z)'];
tmp.Cone = rome_constants.SOC;
rome_constraint(tmp);

% inform the model that a bound has been applied
h.bBoundUsed = true;

% % HACK
% % get the current model
% h = rome_get_current_model();
% 
% % get the primitive uncertain variables for y
% [proj_covar, F_mix] = h.get_cov(y);
% z_std = chol(proj_covar); % projected covariance
% 
% % making bound object
% y_mean = strip_rand(y);
% 
% u = norm2([y_mean'; z_std*strip_certain(y)])';
% bound_obj = 0.5* (y_mean + u);
% 
% % inform the model that a bound has been applied
% h.bBoundUsed = true;
% % END HACK



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
