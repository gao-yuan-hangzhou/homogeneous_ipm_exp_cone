function y = quadexpcone(a, b, mu, k, L)

% QUADEXPCONE Returns y such that mu*exp(a/mu + (b/mu).^2) <= y. 
%
% Currently only accepts scalar a, b, and mu. 
% 
% Uses a Taylor Expansion of exponential function to represent in SOC form.
% k and L control accuracy of approximation. Default Values: k = 4, L = 6.
%
% Note that as part of the approximation, the MATLAB function nchoosek is 
% called and this function will give warnings if k > 26. Also, since this function
% takes 2^L, due to numerical precision issues, L should not be too large.
% Typically, L < 10 should not yield significant numerical errors.
% 
% Calls expcone and hypoquad internally.
%
% Modified By:
% 1. Joel

if(nargin < 5)
    L = 6;
end
if(nargin < 4)
    k = 4;
end

% do scalar expansion
if(~isscalar(a))
    out_sz = size(a);
elseif(~isscalar(b))
    out_sz = size(b);
else
    out_sz = size(mu);
end

% vectorize everyone
a = a(:);
b = b(:);
mu = mu(:);

% get size for vectorization
len = prod(out_sz);

% make new variables
d = rome_model_var(len, 'Cone', rome_constants.NNOC);
x = rome_model_var(len);

% apply constraints
rome_constraint(a + d <= x);
q = hypoquad(mu, d);
rome_constraint(b <= q);

% % use for loop for xm
% S = struct('type', '()', 'subs', []);
% y = [];
% for ii = 1:x.TotalSize
%     S.subs = {ii};
%     y = [y; expcone(x.subsref(S), mu.subsref(S), k, L)];
% end
y = expcone(x, mu, k, L);

% reshape y to mirror input size
y = reshape(y, out_sz);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
