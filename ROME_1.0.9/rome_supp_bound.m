function [bound_obj, S] = rome_supp_bound(y)

% ROME_SUPP_BOUND Creates the support bound for E()^+ in epigraph form
%
% inf_s {max {
%             sup_{z, z_mean} (y_0 + y'z + s'(z_mean - z)),
%             sup_{z, z_mean} (            s'(z_mean - z)) 
%             }}
%
% Input  : y (should be an LDR rome_var)
% Returns: Bounding variable (and optionally s)
%
% Modification History: 
% 1. Joel (11 Feb 2009)

% get the current model
h = rome_get_current_model();

% get the primitive uncertain variables for y
z = h.get_rand_vars();

% make the dual variable S and the bounding variable
S = rome_model_var(y.TotalSize, z.TotalSize);
bound_obj = rome_model_var(y.TotalSize, 'Cone', rome_constants.NNOC); % is nonnegative necessary ??

% apply the bounds
Rhs = bound_obj - S*(z.mean - z);
rome_constraint(y <= Rhs);
rome_constraint(0 <= Rhs);

% inform the model that a bound has been applied
h.bBoundUsed = true;

% % HACK
% % forward bound
% S = rome_model_var(size(y), 'Cone', rome_constants.NNOC);
% y1 = rome_clone(y);
% rome_constraint(S >= y1);
% 
% % backward bound
% R = rome_model_var(size(y), 'Cone', rome_constants.NNOC);
% y2 = rome_clone(y);
% rome_constraint(R >= -y2);
% 
% % convolve
% rome_constraint(y == y1 + y2);
% bound_obj  = S + R + strip_rand(y2);
% % END HACK

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
