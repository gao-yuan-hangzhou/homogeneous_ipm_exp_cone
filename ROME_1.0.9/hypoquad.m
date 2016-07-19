function w = hypoquad(x, y, varargin)

% HYPOQUAD returns a new variable w, such that w^2 <= x*y.
%
% This uses the result that w^2 <= x*y iff 
%
% [w ; (x-y)/2] <=(SOC) (x+y) / 2
%
% OR
%
% [(x+y)/2 ; w ; (x-y)/2] \in SOC
%
% Modified by:
% 1. Joel

% size check
if(isscalar(x))
    out_sz = size(y);
else
    % Should go here if y is a scalar. 
    % If not, or if size is inconsistent, will have an error later
    out_sz = size(x);
end

% vectorize
x = x(:)';
y = y(:)';

% make the new variable if unspecified
if(isempty(varargin))
%     w = rome_model_var([1, prod(out_sz)], 'Cone', rome_constants.NNOC); % is non-negative necessary?
    w = rome_model_var([1, prod(out_sz)]); % is non-negative necessary?
else
    w = varargin{1};
    w = w(:)'; % make into a row vector
end
% w = rome_model_var(size(x)); % is non-negative necessary?

% manually apply SOC constraint (notice implicit vectorization)
q = [(x + y) ./ 2; w; (x - y) ./ 2];
q.Cone = rome_constants.SOC;
rome_constraint(q);

% reshape output back in original shape
w = reshape(w, out_sz);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
