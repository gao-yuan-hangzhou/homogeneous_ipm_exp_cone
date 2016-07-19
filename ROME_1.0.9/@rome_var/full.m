function A = full(A)

% PROF_VAR\FULL Overrides full functionality by converting AffineMap to
% a full matrix. If A is already a full matrix, ignore
%
%   A_full = full(A);
%
%
% Modification History: 
% 1. Joel 
% out = builtin('subsref',A,S);

if(issparse(A))
    A.BiAffineMap = full(A.BiAffineMap);
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
