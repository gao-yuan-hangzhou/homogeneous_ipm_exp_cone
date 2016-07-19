function A = sparse(A)

% PROF_VAR\SPARSE Overrides builtin sparse functionality by converting
% AffineMap of A into a sparse matrix. If A is already a sparse matrix,
% ignore.
%
%   A_sparse = sparse(A);
%
%
% Modification History: 
% 1. Joel 

if(~issparse(A))
    A.BiAffineMap = sparse(A.BiAffineMap);
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
