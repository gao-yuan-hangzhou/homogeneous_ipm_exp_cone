function Y = blk_straighten(X, p, q, rc_flag)
    
% BLK_STRAIGHTEN Straightens a matrix, block-wise
%
% X = [ A B
%       C D ];
%
%   into
%
%       Y = [ A
%             C
%             B
%             D ];      (col)
%
%   or 
%
%       Y = [ A B C D ]; (row)
%
%
% USAGE:
%   Y = blk_straighten(X, p, q)
%   Y = blk_straighten(X, p, q, option) 
%
%   p is the number of columns of each submatrix (A, B, C, D)
%   q is the number of rows in each submatrix (A, B, C, D)
%   option must be either 'row' or 'col'. If omitted, uses 'col' by default
%    
%
% Modified By: 
% 1. Joel
if(nargin == 3)
    rc_flag = 'col';
end

[m, n] = size(X);

if(strcmp(rc_flag, 'col'))
    if(~issparse(X))    
        Y = reshape( X, [ m q n/q ] );
        Y = permute( Y, [ 1 3 2 ] );
        Y = reshape( Y, [ m*n/q q ] );
    else
        col_reps = floor(n ./ q);

        % make permutation index and do block transpose
        perm_ind = reshape(1:n, q, col_reps).';
        Y = reshape(X(:, perm_ind(:)), m * col_reps, q);
    end
elseif(strcmp(rc_flag, 'row'))
    if(~issparse(X))
        Y = reshape( X, [ p m/p n ] );
        Y = permute( Y, [ 1 3 2 ] );
        Y = reshape( Y, [ p m*n/p ] );
    else
        % do block straigten
        Y = blk_straighten(X', q, p, 'col')';
    end    
else
    error('blk_straighten:UnknownOption', 'rc_flag must be either ''row'' or ''col''');
end



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
