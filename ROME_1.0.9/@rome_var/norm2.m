function t = norm2(x)

% ROME_VAR\NORM2 Implements Euclidean Norm for rome_var objects
%
%   v = norm2(x) returns an upper bound on the L2-norm of x. 
% 
%   This function has been 'vectorized' in the following way:
%   If x is a N1 x N2 x N3 ...etc multidimensional array, 
%   the returned value of t will be a 1 x N2 x N3 ... etc multidimensional
%   array such that each element of t provides an upper bound to the
%   corresponding norm of the columns of x as:
% 
%   t(1, n1, n2, n3 ...) >= norm2 (x(:, n1, n2, n3 ....))
%
% Modification History: 
% 1. Joel 

% allow this for vectors only first
% if(~isvector(x))
%     error('rome_var:norm2:InvalidArg', 'Norm2 only accepts vector arguments');
% end

% if x is a row vector, convert into a column vector
sz_x = size(x);
if(isvector(x) && sz_x(1) == 1)
    x = x';
end

% test
if(x.IsRand)
    t = rome_rand_model_var([1, sz_x(2:end)], 'Cone', rome_constants.NNOC);
else
    % make a new model variable
    t = rome_model_var([1, sz_x(2:end)], 'Cone', rome_constants.NNOC);
end

% concatenate the variable
y = [t; x];

% apply SOC restriction
y.Cone = rome_constants.SOC;

% add the constraint
rome_constraint(y);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
