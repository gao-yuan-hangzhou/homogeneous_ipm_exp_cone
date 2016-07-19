function S = outer_sum(x, y)

% OUTER_SUM: Performs pairwise sum of elements of two vectors
%
% USAGE: S = outer_sum(x, y)
%
% x should be a N x 1 column vector
% y should be a 1 x M row vector
%
% Output will be S, an N x M matrix, such that the ij^th element is given
% by: s_ij = x_i + y_j. Notice that the ordering of elements output follows
% the ordering of the first vector (x) supplied.

sz_X = size(x);
sz_Y = size(y);

if(sz_X(2) ~= 1)
    error('outer_sum:FirstColumnVector', 'First Input must be a column vector')
end

if(sz_Y(1) ~= 1)
    error('outer_sum:SecondRowVector', 'Second Input must be a row vector')
end

S = bsxfun(@plus, x, y);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
