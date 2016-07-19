function val = issparse(obj)

% PROF_VAR\ISSPARSE Overrides builtin issparse functionality by checking
% whether the AffineMap of obj is a sparse matrix
%
%   is_sparse_flag = issparse(A);
%
%
% Modification History: 
% 1. Joel 

if(isempty(obj.BiAffineMap))
    val = 1;
else
    val = issparse(obj.BiAffineMap);
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
