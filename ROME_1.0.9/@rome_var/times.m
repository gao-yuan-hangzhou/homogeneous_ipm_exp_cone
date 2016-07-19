function out_var_obj = times(A, B)

% ROME_VAR\TIMES Implements elementwise array multipication in ROME
%
%   C = A .* B 
%
%   At least one of A, B must be a rome_var (cannot be both) 
%   A, B must be of the same size.
%
% Modification History: 
% 1. Joel 

% operation on an empty matrix will yield the other
if(isempty(A))
    out_var_obj = B;
    return;
elseif(isempty(B))
    out_var_obj = A;
    return;
end


% Throw an error if the sizes do not match.
if(~isscalar(A) && ~isscalar(B) && any(size(A) ~= size(B)))
    error('rome_var:times:InvalidArg', 'Input matrices must have the same size.');
end

% Allocate space for output variable
if(isscalar(A))
    out_var_obj = rome_var(size(B));
else
    out_var_obj = rome_var(size(A));
end

% Case 1: A is a (constant) numeric matrix and B is a rome_var
if(isnumeric(A))
    if(isscalar(A))
        out_var_obj.BiAffineMap = spdiags(A * ones(B.TotalSize, 1), 0, B.TotalSize, B.TotalSize) * B.BiAffineMap;
        out_var_obj.DiagMult = A .* B.DiagMult;
    elseif(isscalar(B))
        out_var_obj.BiAffineMap = A(:) * B.BiAffineMap; % computes the outer product
    else
        out_var_obj.BiAffineMap = spdiags(A(:), 0, B.TotalSize, B.TotalSize) * B.BiAffineMap;
    end
    out_var_obj.NumUnmappedVars = B.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = B.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars = B.NumMappedRandVars;
    
% Case 2: A is a rome_var and B is a (constant) numeric matrix
elseif(isnumeric(B))
    if(isscalar(B))
        out_var_obj.BiAffineMap = spdiags(B * ones(A.TotalSize, 1), 0, A.TotalSize, A.TotalSize) * A.BiAffineMap;
        out_var_obj.DiagMult = B .* A.DiagMult;
    elseif(isscalar(A))
        out_var_obj.BiAffineMap = B(:) * A.BiAffineMap; % computes the outer product
    else
        out_var_obj.BiAffineMap = spdiags(B(:), 0, A.TotalSize, A.TotalSize) * A.BiAffineMap;
    end
    out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars = A.NumMappedRandVars;

else
    % Case 3: Mutliplying rome_vars is allowed if 1 is pure certain and 1 is
    % pure uncertain
    if(A.IsCertain && B.IsRand)
        A_expand = repmat(A.BiAffineMap, 1, B.NumMappedRandVars + 1);
        B_expand = reshape(repmat(B.BiAffineMap, A.NumMappedVars + 1, 1), B.TotalSize, size(A_expand, 2));
        
        % handles scalar case by binary singleton expansion
        out_var_obj.BiAffineMap = bsxfun(@times, A_expand, B_expand);
    elseif(A.IsRand && B.IsCertain)
        B_expand = repmat(B.BiAffineMap, 1, A.NumMappedRandVars + 1);
        A_expand = reshape(repmat(A.BiAffineMap, B.NumMappedVars + 1, 1), A.TotalSize, size(B_expand, 2));
        
        % handles scalar case by binary singleton expansion
        out_var_obj.BiAffineMap = bsxfun(@times, A_expand, B_expand);
    else
        % Final Case: Disallowed: Cannot multiply two rome_var objects
        error('rome_var:times:InvalidArg', 'Cannot mutiply two objects of type rome_var');
    end
    
    % Assign outputs
    out_var_obj.NumUnmappedRandVars = max(A.NumUnmappedRandVars, B.NumUnmappedRandVars);
    out_var_obj.NumUnmappedVars = max(A.NumUnmappedVars, B.NumUnmappedVars);
    out_var_obj.NumMappedRandVars = max(A.NumMappedRandVars, B.NumMappedRandVars);
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
