function obj = rome_empty_var(varargin)

% ROME_EMPTY_VAR: Makes an empty rome_var object, can be used as a
% pre-allocator

h = rome_get_current_model();

% only accept numeric (size) argument
if(any(~cellfun(@isnumeric, varargin)))
    error('rome_empty_var:OnlyNumericArgs', 'rome_empty_var only accepts numeric arguments');
end
obj = rome_var(varargin{:});

try
    obj.NumUnmappedVars = h.NumVars;
catch ME
    % Will Enter here if the user forgets to call rome_model beforehand
    if(strcmp(ME.identifier, 'MATLAB:nonStrucReference'))
        warning('rome_empty_var:ctor:UninitializedModel', ...
            'No Rome Model found: Creating a new model. To remove this error message, call rome_model first');
        h = rome_model;
    else
        rethrow(ME);
    end
end
obj.BiAffineMap = spalloc(obj.TotalSize, 1, 0); 
obj.NumUnmappedRandVars = h.NumRandVars;

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
