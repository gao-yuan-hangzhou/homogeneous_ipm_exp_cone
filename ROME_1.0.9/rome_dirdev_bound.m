function bound_obj = rome_dirdev_bound(y)

% ROME_DIRDEV_BOUND Creates the directional deviation bound for E()^+ in epigraph form
%
% Input  : y (should be an LDR rome_var)
% Returns: Bounding variable
%
% Modification History: 
% 1. Joel (13 Feb 2009)

% get the current model
h = rome_get_current_model();

% vectorize
orig_sz = size(y);
y = y(:);
len = length(y);

% zhat = h.get_rand_vars(y.mean);
% z = h.get_rand_vars(y);
z = h.get_rand_vars();
zhat = z.mean;

% get directional deviations for y
proj_fdev = h.get_dirdev(y, 1);
[proj_bdev, affine_sigma_mix] = h.get_dirdev(y, -1);

proj_fdev = repmat(proj_fdev', 2*len, 1);
proj_bdev = repmat(proj_bdev', 2*len, 1);

F = affine_sigma_mix(:, 2:end);
g = affine_sigma_mix(:, 1);

% don't use linear rules - do this directly: faster
y0   = strip_rand(y);
y_vec = strip_certain(y);

% make aux variable
s0 = rome_model_var(len);
s_vec  = rome_model_var(len, size(F, 1));

% make a new aux variable
x0 = rome_model_var(len);
x_vec  = rome_model_var(len, size(F, 1));

% make aux t
t = rome_model_var(len);

% apply constraints
rome_constraint(y.mean - (s0 + s_vec * (F * zhat + g)) <= t);
rome_constraint(x0 + x_vec * g == y0); 
rome_constraint(x_vec * F == y_vec'); 

% make aux variables
u = rome_model_var(2*len, size(F, 1));

% rome_constraint(u1 >=  proj_fdev .* (s_vec - x_vec));
% rome_constraint(u1 >= -proj_bdev .* (s_vec - x_vec));
% rome_constraint(u2 >=  proj_fdev .* s_vec);
% rome_constraint(u2 >= -proj_bdev .* s_vec);

rome_constraint(u >=  proj_fdev .* [s_vec - x_vec; s_vec]);
rome_constraint(u >= -proj_bdev .* [s_vec - x_vec; s_vec]);

% create the conic bound
lambda1 = rome_model_var(len, 'Cone', rome_constants.NNOC);
lambda2 = rome_model_var(len, 'Cone', rome_constants.NNOC);

u1 = u(1:len, :);
u2 = u(len+(1:len), :);

% make psi function
z_sigmahat = F*zhat + g;
b1 = rome_model_var(len);
b2 = rome_model_var(len);
rome_constraint(s0 - x0 + (s_vec - x_vec)*z_sigmahat <= b1);
rome_constraint(s0 + s_vec * z_sigmahat <= b2);

q1 = quadexpcone(b1, 1/sqrt(2) * norm2(u1')', lambda1);
q2 = quadexpcone(b2, 1/sqrt(2) * norm2(u2')', lambda2);

% make the output variable
bound_obj = reshape(t + exp(-1) * (q1 + q2), orig_sz); 

% inform the model that a bound has been applied
h.bBoundUsed = true;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
