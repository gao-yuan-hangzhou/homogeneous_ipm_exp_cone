function [zz Map] = rome_clone_uncertainty(z)
% ROME_CLONE_UNCERTAINTY Clones a random variable for specifying
% conditional constraints
%
% (A) Variables are specified by a rome_var
% [zz Map] = rome_clone_uncertainty(z)
%
% (B) Variables are specified by an index set (used in internal calls)
% [zz Map] = rome_clone_uncertainty(ind)
%
% Input  :
%   In case (A), z must be a random and primitive rome_var.
%   In case (B), ind must be a set of indices to random model variables.
% Returns:
%   zz is the cloned random variable with the same size as z.
%   Map is a struct representing the mapping relationship between z and zz.


% case 1: input is a rome_var
if (isa(z, 'rome_var'))
    % check if z is random
    if  (~z.IsRand)
        error('rome_clone_uncertainty:InvalidInput', ...
            'Input is not purely random. Input must be either a purely random and primitive rome_var, or an index set.')
    end
    
    % check if z is primitive
    zMap = logical(z.BiAffineMap(:, 2:end)); % ignoring the constant term
    if any(sum(zMap, 2) ~= 1)
        error('rome_clone_uncertainty:InvalidInput', ...
        'Input is not primitive. Input must be either a purely random and primitive rome_var, or an index set.')
    end
    
    % find the primitive index set
    [c r] = find(zMap);
    ind = zeros([1 z.TotalSize]); % the index set is a row vector
    ind(c) = r + z.NumUnmappedRandVars;
    
    % create a clone of z
    zz = rome_rand_model_var(size(z));
    
% case 2: input is an index set
elseif isnumeric(z)
    % force the index set to be a row vector 
    ind = z(:)';
    % create a clone rand var as a column vector
    zz = rome_rand_model_var([numel(ind) 1]);
    
else
    error('rome_clone_uncertainty:InvalidInput', ...
        'Input must be either a purely random and primitive rome_var, or an index set.')
end

% record the mapping relationship (a shift in random variable index)
Map.origPrimInd = ind;
Map.clonedPrimInd = zz.NumUnmappedRandVars + (1:zz.NumMappedRandVars);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.