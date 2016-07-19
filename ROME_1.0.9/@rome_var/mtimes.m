function out_var_obj = mtimes(A, B)

% ROME_VAR\MTIMES Implements matrix multiplication in ROME
%
%   C = A * B 
%
%   At least one of A, B must be a rome_var
%
% Modification History: 
% 1. Joel 

if(isscalar(A) || isscalar(B))
    out_var_obj = A .* B;
    return;
else
    % Error check dimensions
    if(size(A, 2) ~= size(B, 1))
        error('rome_var:mtimes:InvalidArg', 'Inner matrix dimensions must agree.');
    end

    % Notice that size(B) has been overloaded.
    out_var_obj = rome_var(size(A, 1), size(B, 2));
end

% Case 1: A is a (constant) numeric matrix and B is a rome_var
if(isnumeric(A))
%     if(isscalar(A) || isscalar(B))
%         out_var_obj.BiAffineMap = A(:) * B.BiAffineMap;
%     else
        num_cols = size(B.BiAffineMap, 2); % define number of columns of B
        wrapped_B = reshape(B.BiAffineMap, size(A, 2), size(B, 2) * num_cols);  % wrap
        out_var_obj.BiAffineMap = reshape(A * wrapped_B, out_var_obj.TotalSize, num_cols); % unwrap
%     end
    out_var_obj.NumUnmappedVars = B.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = B.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars = B.NumMappedRandVars;

% Case 2: A is a rome_var and B is a (constant) numeric matrix
elseif(isnumeric(B))
%     if(isscalar(B) || isscalar(A))
%         out_var_obj.BiAffineMap = B(:) * A.BiAffineMap;
%     else
        num_cols = size(A.BiAffineMap, 2);  % Define number of columns of A
        wrapped_A = reshape(A.BiAffineMap', size(A, 1)* num_cols, size(B, 1)); % transposition on A
        out_var_obj.BiAffineMap = reshape(wrapped_A * B, num_cols, out_var_obj.TotalSize)'; % notice final transposition here
%     end
    out_var_obj.NumUnmappedVars = A.NumUnmappedVars;
    out_var_obj.NumUnmappedRandVars = A.NumUnmappedRandVars;
    out_var_obj.NumMappedRandVars = A.NumMappedRandVars;

else
    % Case 3: Mutliplying rome_vars is allowed if 1 is pure certain and 1 is
    % pure uncertain
    if(A.IsCertain && B.IsRand)
        % define some variables for convenience
        p = size(A, 1);  % num rows of A (= num output rows)
        m = size(B, 1);  % inner size (= num cols of A = num rows of B)
        q = size(B, 2);  % num cols of B (= num output cols)
        N = A.NumMappedVars + 1;     % Internal Num Cols of A
        M = B.NumMappedRandVars + 1; % Internal Num Cols of B
        
        % perm ind for A
        A_perm_ind = reshape(1:A.TotalSize, p, m).';
        A_perm = A.BiAffineMap(A_perm_ind(:), :);     % permute rows of A
        A_perm = blk_transpose(A_perm, m, N) .';      % block transpose followed by transpose of A

        % preprocess B
        B_perm = blk_transpose(B.BiAffineMap, m, M);
        
        % output
        out_perm = blk_straighten(A_perm * B_perm, N, M);
        out_perm = blk_transpose(out_perm, N, M, 'inner');  % do a final inner transpose
        out_var_obj.BiAffineMap = reshape(out_perm.', N * M, p*q).';
    elseif(A.IsRand && B.IsCertain)
        % define some variables for convenience
        p = size(A, 1);  % num rows of A (= num output rows)
        m = size(B, 1);  % inner size (= num cols of A = num rows of B)
        q = size(B, 2);  % num cols of B (= num output cols)
        N = A.NumMappedRandVars + 1; % Internal Num Cols of A
        M = B.NumMappedVars + 1; % Internal Num Cols of B
        
        % perm ind for A
        A_perm_ind = reshape(1:A.TotalSize, p, m).';
        A_perm = A.BiAffineMap(A_perm_ind(:), :);     % permute rows of A
        A_perm = blk_transpose(A_perm, m, N) .';      % block transpose followed by transpose of A

        % preprocess B
        B_perm = blk_transpose(B.BiAffineMap, m, M);
        
        % output
        out_perm = blk_straighten(A_perm * B_perm, N, M);
        out_var_obj.BiAffineMap = reshape(out_perm.', N * M, p*q).';
    else
        % Final Case: Disallowed: Cannot multiply two rome_var objects
        error('rome_var:mtimes:InvalidArg', 'Cannot mutiply two objects of type rome_var');
    end
    
    % Assign outputs
    out_var_obj.NumUnmappedRandVars = max(A.NumUnmappedRandVars, B.NumUnmappedRandVars);
    out_var_obj.NumUnmappedVars = max(A.NumUnmappedVars, B.NumUnmappedVars);
    out_var_obj.NumMappedRandVars = max(A.NumMappedRandVars, B.NumMappedRandVars);
    
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
