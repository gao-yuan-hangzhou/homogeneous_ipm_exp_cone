function out_var_obj = ctranspose(A)

% ROME_VAR\CTRANSPOSE Implements complex conjugate transpose in ROME
%
%   C = A'
%
%   Same as real transposition
%
% Modification History: 
% 1. Joel 

% Throw an error if not a matrix
if(length(size(A)) > 2)
    error('rome_var:ctranspose:InvalidArg', 'Ctranspose only operates on 2 dimensional matrices');
end

% Allocate space for output variable
out_var_obj = rome_var(size(A, 2), size(A, 1));

% create permutation matrix
perm_ind = reshape(1:A.TotalSize, size(A, 1), size(A, 2))';
perm_ind = perm_ind(:);

perm_mat = speye(A.TotalSize);
perm_mat = perm_mat(perm_ind, :);

% mutiply
out_var_obj.BiAffineMap = perm_mat * A.BiAffineMap;

% Define output
out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
out_var_obj.NumMappedRandVars = A.NumMappedRandVars;

% transposition generally destroys diagonality, unless it's a vector
if(isvector(A))
    out_var_obj.DiagMult = A.DiagMult;
end
out_var_obj.Cone = A.Cone;
out_var_obj.Continuity = A.Continuity;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
