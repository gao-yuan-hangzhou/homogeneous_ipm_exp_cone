function q = norm1(x)

% ROME_VAR\NORM1 Implements Manhattan Norm for rome_var objects
%
%   v = norm1(x) returns an upper bound on the L1-norm of x. 
% 
%   This function has been 'vectorized' in the folloing way:
%   If x is a N1 x N2 x N3 ...etc multidimensional array, 
%   the returned value of t will be a 1 x N2 x N3 ... etc multidimensional
%   array such that each element of t provides an upper bound to the
%   corresponding norm of the columns of x as:
% 
%   q(1, n1, n2, n3 ...) >= norm2 (x(:, n1, n2, n3 ....))
%
% Modification History: 
% 1. Joel 

if(x.IsCertain)
    % make a new model variable
    t = rome_model_var(size(x), 'Cone', rome_constants.NNOC);
elseif(x.IsRand)
    % make a new model rand variable
    t = rome_rand_model_var(size(x), 'Cone', rome_constants.NNOC);
else
    error('rome_var:norm1:NoLDR', 'Norm1 does not accept LDR inputs');
end

% apply constraints
rome_constraint(x <= t);
rome_constraint(-t <= x);

% returns the sum 
q = sum(t, 1);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
