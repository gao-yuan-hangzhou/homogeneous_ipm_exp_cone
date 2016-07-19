function [bound_obj, S] = rome_pwl_supp_bound(y, b_vec, a_vec)

% ROME_PWL_SUPP_BOUND Generates the support bound for a piece-wise linear 
% function of an LDR. (TESTING)
%
% inf_s {
%   max_k {
%       sup_{z, z_mean} (b_k + a_k y_0 + a_k y'z + s'(z_mean - z))
% }}
%
% coeff should be formatted as 
%   coeff = [b_1, a_1; b_2, a_2; b_3, a_3 ... ; b_k, a_k] 
%             
% 
% Modification History: 
% 1. Joel (7 Sept 2009)

% get the current model
h = rome_get_current_model();

% get the primitive uncertain variables for y
z = h.get_rand_vars();

% make the dual variable S and the bounding variable
S = rome_model_var(y.TotalSize, z.TotalSize);
bound_obj = rome_model_var(y.TotalSize, 'Cone', rome_constants.NNOC); % is nonnegative necessary ??

% apply the bounds
Rhs = bound_obj - S*(z.mean - z);
for ii = 1:size(b_vec, 2)
    if(ndims(a_vec) > 2)
        rome_constraint(b_vec(:, ii) + a_vec(:, :, ii)* y <= Rhs);
    else
        % optimization to handle case of elements of a_vec being scalar
        rome_constraint(b_vec(:, ii) + a_vec(ii) * y <= Rhs);
    end
end

% inform the model that a bound has been applied
h.bBoundUsed = true;



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
