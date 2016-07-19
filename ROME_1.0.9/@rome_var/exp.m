function y = exp(x, k, L)

% ROME_VAR\EXP Returns y such that exp(x) <= y. 
%
% Currently only accepts scalar x and y. 
% Uses a Taylor Expansion of exponential function to represent in SOC form. 
%
%   exp(x) = c + (sum{i = 1:k} (a_i * (x/2^L) + b_i))^2^L
% 
% Default Values: k = 4, L = 6.
%
% Note that as part of the approximation, the MATLAB function nchoosek is 
% and this function will give warnings if k > 26. Also, since this function
% takes 2^L, due to numerical precision issues, L should not be too large.
% Typically, L < 10 should not yield significant numerical errors.
% 
%
% Modified By:
% 1. Joel

if(~isscalar(x))
    error('rome_var:exp:ScalarArgsOnly', 'expcone only accepts scalar x');
end

% default vals
if(nargin < 3)
    L = 6;
end
if(nargin < 2)
    k = 4;
end

y = expcone(x, 1, k, L);

% % call a helper function to obtain exponential cone coefficients
% [a, b, c] = expconecoef(k);
% 
% % perform expansion
% F = 2^L;
% z = (a / F) * x + b;    % overloading power
% 
% % structure for subscripted referencing
% S = struct('type', '()', 'subs', []);
% 
% sum_z = c;
% for ii = 1:k
%     S.subs = {ii};
%     sum_z = sum_z + z.subsref(S) .^ (2*ii);
% end
% 
% y = sum_z^F;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
