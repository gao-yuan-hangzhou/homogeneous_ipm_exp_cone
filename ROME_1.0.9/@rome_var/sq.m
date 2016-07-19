function y = sq(x, Q)

% ROME_VAR\SQ Implements "Square" Function, y = x' * Q * x
%
% This function is implemented to implement and optimize squaring. Notice
% that Q supplied should be positive semidefinite
%
%   y = sq(x);
%   y = sq(x, Q);
%
%  ||Q ^ (1/2) x||^2 <= t can be expressed as: sq(x, Q) <= t
%
%
% Modification History: 
% 1. Joel 

% allow this for vectors only first
if(~isvector(x))
    error('rome_var:sq:VectorArgOnly', 'SQ only accepts vector arguments');
end

if(nargin < 2)
    Q = speye(length(x));
end

% make a new model variable
y = rome_model_var('Cone', rome_constants.NNOC);

% WORKING
% concatenate the variable
z = [(y+1) ./ 2; x(:); (y-1) ./ 2];
% END WORKING

% % TESTING
% % concatenate the variable
% z = [y; x(:)];
% 
% % apply Non-linearity
% z.NonLinearity = 'SQ';
% z.ExtraData    = Q; 
% % END TESTING

% apply SOC restriction
z.Cone = rome_constants.SOC;

% add the constraint
rome_constraint(z);



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
