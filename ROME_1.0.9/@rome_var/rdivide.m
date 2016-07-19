function out_var_obj = rdivide (A, B)

% ROME_VAR\TIMES Implements elementwise array (right) standard division in ROME
%
%   C = A ./ B 
%
%   A must be a rome_var
%   B cannot be a rome_var
%   A, B must be of the same size.
%
% Modification History: 
% 1. Joel 

% operation on an empty matrix will yield the other
if(isempty(A) || isempty(B))
    error('rome_var:rdivide:InvalidArg', 'Input matrices cannot be empty.');
end


% Throw an error if the sizes do not match.
if(~isscalar(A) && ~isscalar(B) && all(size(A) ~= size(B)))
    error('rome_var:rdivide:InvalidArg', 'Input matrices must have the same size.');
end

% Allocate space for output variable
if(isscalar(A))
    out_var_obj = rome_var(size(B));
else
    out_var_obj = rome_var(size(A));
end

% Case 1: A is a (constant) double matrix and B is a rome_var
if(isnumeric(A))
    error('rome_var:rdivide:InvalidArg', 'B cannot be a rome_constants.');
    
% Case 2: A is a rome_var and B is a (constant) double matrix
elseif(isnumeric(B))
    if(isscalar(B))
        out_var_obj.BiAffineMap = spdiags(ones(A.TotalSize, 1) ./ B(:), 0, A.TotalSize, A.TotalSize) * A.BiAffineMap;
        out_var_obj.DiagMult  = A.DiagMult ./ B;
    elseif(isscalar(A))
        out_var_obj.BiAffineMap = (1./ B(:)) * A.BiAffineMap; % compute outer produc
    else
        out_var_obj.BiAffineMap = spdiags(1./B(:), 0, A.TotalSize, A.TotalSize) * A.BiAffineMap;
    end    
    out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars = A.NumMappedRandVars;
end



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
