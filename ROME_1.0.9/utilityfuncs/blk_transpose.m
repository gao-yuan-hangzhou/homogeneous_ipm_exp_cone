function Y = blk_transpose(X, p, q, inner_flag)
    
% BLK_TRANSPOSE Transposes a matrix, block-wise. Does REAL Transposition
%
% Outer Transposition: Keeps submatrices untransposed. 
%   [A, B   -->  [A, C
%    C, D]  -->   B, D];
%
% Inner Transposition: Leaves blocks and transposes submatrices
%   [A, B   -->  [A.', B.'
%    C, D]  -->   C.', D.'];

% USAGE:
%   Y = blk_transpose(X, p, q)
%   Y = blk_transpose(X, p, q, option) 
%
%   p is the number of rows of each submatrix (A, B, C, D)
%   q is the number of columns in each submatrix (A, B, C, D)
%   option must be either 'inner' or 'outer'. If omitted, uses 'outer' by default
%
%
% Modified By: 
% Joel

% check number of args
if(nargin == 3)
    inner_flag = 'outer';
end

% define size of X
[m, n] = size(X);

% OUTER Transpose
if(strcmp(inner_flag, 'outer'))
    
    if(~issparse(X))
        Y = reshape( X, [ p m/p q n/q ] );
        Y = permute( Y, [ 1 4 3 2 ] );
        Y = reshape( Y, [ p*n/q q*m/p] );
    else
        row_reps = floor(m ./ p);
        col_reps = floor(n ./ q);

        % make permutation index
        perm_ind = reshape(1:n, q, col_reps).';
        Y = X(:, perm_ind(:));
        
        % straighten first
        Y = blk_straighten(Y.', q, p).';
        
        % reshape
        Y = reshape(Y, p * col_reps, q * row_reps);
    end

% INNER transpose
elseif(strcmp(inner_flag, 'inner'))
    if(~issparse(X))
        Y = reshape( X, [ p m/p q n/q ] );
        Y = permute( Y, [ 3 2 1 4 ] );
        Y = reshape( Y, [ q*m/p p*n/q ] );
    else
        Y = blk_transpose(X.', q, p, 'outer');
    end 
else
    error('blk_transpose:InvalidOption', 'Last Argument must be either "inner" or "outer".');
end



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
