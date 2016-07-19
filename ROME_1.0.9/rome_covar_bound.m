function [bound_obj, u] = rome_covar_bound(y)

% ROME_COVAR_BOUND Creates the covar bound for E()^+ in epigraph form
%
% sup_{z_mean}{0.5 (y_0 + y'z_mean) 
%             + 0.5 * sqrt( (y_0 + y'z_mean)^2 + y' * Covar * y ) }
%
% Input  : y (should be an LDR rome_var)
% Returns: Bounding variable
%
% Modification History: 
% 1. Joel (13 Feb 2009)

% get the current model
h = rome_get_current_model();

% get the primitive uncertain variables for y
[proj_covar, F_mix] = h.get_cov(y);
z_std = chol(proj_covar); % projected covariance
%z_std = full(proj_covar)^(0.5); % projected covariance % HACK!
z_std = real(z_std); % for numerical niceness

% make the auxilliary variable u (for the mean)
u = rome_model_var(y.TotalSize);

% bound the mean of y
rome_constraint(y.mean <= u);

% make new variable
X = rome_model_var(size(F_mix, 1), y.TotalSize);
rome_constraint(strip_certain(y) == F_mix' * X); 

% temporary hack (DO NOT COMMIT)
% z_std = sparse(z_std);

% making bound object
bound_obj = 0.5* (u + norm2([u'; z_std * X])');

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
