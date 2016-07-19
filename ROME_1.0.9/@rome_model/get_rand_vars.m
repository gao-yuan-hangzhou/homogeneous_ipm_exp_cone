function obj = get_rand_vars(model_obj, ind)

% ROME_MODEL\GET_VARS Returns a (uncertain) rome_var which 
% indexes into the global uncertain variable space based on ind. 
%
% ind should be a vector of indices into the global vector of variables
% and best if sorted. Returns a column vector of uncertain rome vars
%
% Modification History: 
% 1. Joel 

% Case 0: no ind specified, return vector of uncertainties (no means)
if(nargin == 1)
   obj = rome_var(model_obj.NumRandVars);
   obj.BiAffineMap = [zeros(model_obj.NumRandVars, 1), speye(model_obj.NumRandVars)];
   obj.DiagMult = 1;
   obj.NumMappedRandVars = model_obj.NumRandVars;
   
   obj = obj(find(~model_obj.ZIsMean)); % don't return means
   return;
end

% Case 1: if ind is a rome_var
if(isa(ind, 'rome_var'))
    if(ind.IsConst || ind.IsCertain)
        % Certain or constant return empty matrix
        obj = [];
        return;
    else
        % LDR or uncertain do a recursive call
        numeric_ind = ind.NumUnmappedRandVars + (1:ind.NumMappedRandVars);
        obj = get_rand_vars(model_obj, numeric_ind);
        obj.DiagMult = 1;
        return;
    end
end

% Case 2: if ind is char (must be 'all') flag
if(ischar(ind))
    if(strcmpi(ind, 'all'))
        obj = rome_var(model_obj.NumRandVars);
        obj.BiAffineMap = [zeros(model_obj.NumRandVars, 1), speye(model_obj.NumRandVars)];
        obj.DiagMult = 1;
        obj.NumMappedRandVars = model_obj.NumRandVars;
    else
        error('rome_model:get_rand_vars:IncorrectInd', ...
            'Second argument must be the string ''all''.');
    end
end

% Case 3: if ind is numeric
% throw an exception if ind is out of bounds
min_ind = min(ind(:));
max_ind = max(ind(:));

if((min_ind <= 0) || (max_ind > model_obj.NumRandVars))
    error('rome_model:get_rand_vars:IndexOutOfBounds', ...
        'Indices must be positive numbers no greater than the number of uncertain variables in the model');
end

% make the object
obj = rome_var(numel(ind));
obj.NumUnmappedRandVars = min_ind - 1;
obj.NumMappedRandVars = max_ind - min_ind + 1;

% use spconvert to find AffineMap
row_ind = (1:numel(ind))';
col_ind = 1 + ind(:) - obj.NumUnmappedRandVars;

% make the biaffine map
obj.BiAffineMap = spconvert([row_ind, col_ind, ones(numel(ind), 1)]);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
