function obj = get_vars(model_obj, ind)

% ROME_MODEL\GET_VARS Returns a (deterministic) rome_var which 
% indexes into the global variable space based on ind. 
%
% ind should be a vector of indices into the global vector of variables
% and best if sorted. Returns a column vector of prov_
%
% Modification History: 
% 1. Joel 

% return all if no ind specified
if(nargin == 1)
   obj = rome_var(model_obj.NumVars);
   obj.BiAffineMap = [zeros(model_obj.NumVars, 1), speye(model_obj.NumVars)];
   obj.DiagMult = 1;
   return;
end

% throw an exception if ind is out of bounds
min_ind = min(ind(:));
max_ind = max(ind(:));

if((min_ind <= 0) || (max_ind > model_obj.NumVars))
    error('rome_model:get_vars:IndexOutOfBounds', ...
        'Indices must be positive numbers no greater than the number of variables in the model');
end

% make the object
obj = rome_var(numel(ind));
obj.NumUnmappedVars = min_ind - 1;

% use spconvert to find AffineMap
row_ind = (1:numel(ind))';
col_ind = 1 + ind(:) - obj.NumUnmappedVars;

% make the biaffine map
obj.BiAffineMap = spconvert([row_ind, col_ind, ones(numel(ind), 1)]);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
