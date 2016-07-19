function y = expcone(x, mu, k, L)

% EXPCONE Returns y such that mu*exp(x/mu) <= y. 
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

% size check
if(isscalar(x))
    out_sz = size(mu);
else
    % should go here if mu is a scalar. 
    % If not, or if size is inconsistent, will have an error later
    out_sz = size(x);
end

% if(~isscalar(x) || ~isscalar(mu))
%     error('rome_var:expcone:ScalarArgsOnly', 'expcone only accepts scalar x and mu');
% end

% default vals
if(nargin < 4)
    L = 6;
end
if(nargin < 3)
    k = 4;
end

% vectorize
x = x(:)';
mu = mu(:)';
len = prod(out_sz);

% OLD CODE 
% % make output variable
% y = rome_model_var(out_sz, 'Cone', rome_constants.NNOC);

% reduction on powers of L
% rhs = y;
% for ii = 1:L
%     rhs = hypoquad(mu, rhs);
% end

% temp (optimized reduction on powers of L)
w = rome_model_var([L+1, len], 'Cone', rome_constants.NNOC);
hypoquad(ones(L, 1) * mu, w(1:end-1, :), w(2:end, :));

% call a helper function to obtain exponential cone coefficients
[a, b, c] = expconecoef(k);

% create auxilliary variables (total = k), using overloaded operators for
% scalar rome_vars
F = 2^L;
z = (a / F) * x + b * mu;

% perform expansion of exp function
% sum_z = c * mu ;
% S = struct('type', '()', 'subs', []);
% for ii = 1:k
%     S.subs = {ii};
%     sum_z = sum_z + powercone(z.subsref(S), mu, 2*ii);
% end

% % TEST 1
% for ii = 1:k
%     sum_z = sum_z + powercone(z(ii), mu, 2*ii);
% end

% TEST 2
pow_mat = repmat(2*(1:k)', 1, len); % matrix of powers
sum_z = c * mu + sum(powercone(z, ones(k,1) * mu, pow_mat));

% final bound
rome_constraint(sum_z <= w(end, :));
% rome_constraint(sum_z == rhs);

% assign output
y = reshape(w(1, :), out_sz);

function [a,b,c] = expconecoef(k)

% EXPCONECOEF   Helper function to create expansion coefficients for exp cone 
% in even powers of x
%
%    k = Talor expension of to x^2k
%    Return cone coefficients of the form   c + sum_{k=1:k} (a(k) x + b(k))^(2k)
%    Example:
%       L=5;H = 2^6;
%       [a b c] = expconecoef(L);
%       x=50; 
%       (sum((x/H*a + b).^([2:2:2*L]'))+c)^H   % approximation of exp(x)
%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification History
% Date         Author    Comments
% ----         ------    ----------------------------------------------------
% 22/06/2005   Melvyn    Initial Creation
% 08/10/2008   Joel      Vectorization and comments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 2*k;            % order of original polynomial
a = zeros(k,1);     % output (coefficient on each affine term)
b = zeros(k,1);     % output (constant in each affine term)

% start with original coefficients on exp(x)
tmpcoef = 1./ factorial(0:n);

% iterate 
for ii=k:-1:1
    L = 2*ii;
    a(ii) = tmpcoef(L+1)^(1/L);          % compute coefficent values
    b(ii) = tmpcoef(L)/L/(a(ii)^(L-1));  
    
    % update tmp coef
    tmpcoef(1:L+1) = tmpcoef(1:L+1) - vecbino(L) .* (a(ii).^(0:L)) .* (b(ii).^(L:-1:0));
end

% final constant
c = tmpcoef(1);

%-------------------------
function vec = vecbino(n)
% create binomial coeffieints of n choose k, k=0 to n
vec = zeros(1, n+1);
for ii=0:n
    vec(ii+1)= nchoosek(n,ii);
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
