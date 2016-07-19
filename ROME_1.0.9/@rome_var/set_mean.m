function primitive_mean_val  = set_mean(obj, mean_val)

% ROME_VAR\SET_MEAN Special function to set the numeric mean for a
% uncertain variable. Use this instead of the (less efficient but more general)
% construction:
%
%   rome_constraint(z.mean == mean_val)
%
% when mean_val is a numeric quantity, and z is full-dimensional.
%
% USAGE:
% z = rome_rand_model_var(10, 20);
% z.set_mean(z, ones(10, 20));
%
%
% USE THIS WITH CAUTION. It is your responsibility to verify that z is
% full-dimensional
%
%
% Modification History: 
% 1. Joel (30 Apr 2009)

% check that mean_val is numeric
if(~isnumeric(mean_val))
    error('rome_set_mean:NeedNumeric', 'Supplied Mean value must be numeric.');
end

% check that obj is pure rand
if(~obj.IsRand)
    error('rome_set_mean:NeedRand', 'Supplied Rome Var must be pure uncertain.');
end

% get handle to current model
h_curr_model = rome_get_current_model();

% calculate primitive mean value
b = mean_val - obj.BiAffineMap(:, 1);
A = obj.BiAffineMap(:, 2:end);
primitive_mean_val = A \ b;

% communicate mean value to the model
h_curr_model.set_mean(obj.NumUnmappedRandVars, primitive_mean_val);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
