 function out_mat = insert(obj, in_data, arg)

% ROME_VAR\INSERT Substitutes input values into object
%
%   C1 = A(x_input)             % substitutes into deterministic portion
%   C2 = A(x_input, 'rand')     % substitutes into uncertain portion
%   C3 = A(x_input, z_input);   % substitutes into deterministic (x) and
%                                 uncertain (z) portions
%
% Modification History: 
% 1. Joel 
% out = builtin('subsref',A,S);
if(nargin == 2)
    arg = [];
end
 
if(isempty(obj.BiAffineMap))
    error('rome_var:insert:EmptyArg', 'Input rome_var object must not have an empty map');
end

% if(~isnumeric(in_data))
%     error('rome_var:insert:InvalidArgType', 'Input data must be a numeric type matrix');
% end

if(~isscalar(in_data) && (obj.IsCertain) && (numel(in_data) ~= obj.NumMappedVars))
    error('rome_var:insert:InvalidCertainArg', 'Input certain data must either be a real scalar or match the objects input size');
end

% don't use numel here because we allow insertion of uncertain variables
if(~isscalar(in_data) && (obj.IsRand) && (prod(size(in_data)) ~= obj.NumMappedRandVars))
    error('rome_var:insert:InvalidRandArg', 'Input uncertain data must either be a real scalar or match the objects input size');
end

% perform substitution
% if(obj.IsDiag)
    % Optimization for Diagonal variables
%     out_mat = obj.DiagMult .* in_data(:) + obj.BiAffineMap(:, 1);
%     out_mat =  reshape(out_mat, size(obj));
if(obj.IsCertain || obj.IsRand)
    % Simple case: object is either a pure certain variable or a pure
    % uncertain variable.
    out_mat = obj.BiAffineMap(:, 1) + obj.BiAffineMap(:, 2:end) * in_data(:);
    out_mat =  reshape(out_mat, size(obj));
else
    % More complex case: object is a full LDR
    % Check whether we are only specifying uncertain data
    if(strcmp(arg, 'rand'))
        % reshape affinemap and multiply with uncertain data
        coeff_matrix = reshape(obj.BiAffineMap, obj.TotalSize * (1 + obj.NumMappedVars), 1 + obj.NumMappedRandVars);
        coeff_matrix = coeff_matrix * [1; in_data(:)];
        
        % shape it back to original form
        out_mat = obj;
        out_mat.BiAffineMap = reshape(coeff_matrix, obj.TotalSize, 1+obj.NumMappedVars);
        out_mat.NumMappedRandVars = 0;
        out_mat.NumUnmappedRandVars = 0;
    else
        certain_data  = in_data;
        rand_data = arg;
        
        % when both rand data and certain data are specified
        num_cols = size(obj.BiAffineMap, 2);
        const_ind = 1:(obj.NumMappedVars+1):num_cols;   % index of constant columns
        ind = setdiff(1:num_cols, const_ind);

        int_mat = blk_straighten(obj.BiAffineMap(:, ind), obj.TotalSize, obj.NumMappedVars) * certain_data(:);
        int_mat = obj.BiAffineMap(:, const_ind) + reshape(int_mat, obj.TotalSize, obj.NumMappedRandVars + 1);
        if(~isempty(rand_data))
            out_mat = reshape(int_mat(:, 1) + int_mat(:, 2:end) * rand_data(:), size(obj));
        else
            out_mat = reshape(int_mat, [size(obj), obj.NumMappedRandVars + 1]);
        end
    end
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
