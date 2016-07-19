function obj = rome_model_var(varargin)

% ROME_MODEL_VAR Describes a model level variable in ROME
%
%
% Modification History:
% 1. Joel
%
% check that none of the arguments are rome_var objects. 
if(any(cellfun(@(x) isa(x, 'rome_var'), varargin)))
    error('rome_model_var:NoVarInput', ...
        'rome_model_var does not accept rome_var inputs. Did you mean to call rome_linearrule instead?');
end

% construct the superclass
obj = rome_var(varargin{:});

% Get the current model and query its state
h_curr_model = rome_get_current_model();

try
    % define number of unmapped variables
    obj.NumUnmappedVars = h_curr_model.NumVars;
catch ME
    % Will Enter here if the user forgets to call rome_model beforehand
    if(strcmp(ME.identifier, 'MATLAB:nonStrucReference'))
        warning('rome_model_var:ctor:UninitializedModel', ...
            'No Rome Model found: Creating a new model. To remove this error message, call rome_model first');
        h_curr_model = rome_model;
    else
        rethrow(ME);
    end
end

% make Uncertain variables
obj.NumUnmappedRandVars = 0;
obj.NumMappedRandVars = 0;

% create the mapping rules
obj.BiAffineMap = [spalloc(obj.TotalSize, 1, 1), speye(obj.TotalSize)];

% All model variables are have identity mapping by default
obj.DiagMult = 1;

% add this variable to the current model
h_curr_model.add_new_var(obj);

% Apply the conic constraint at this time. Check for
% presence/absence of constraint done internally in
% rome_constraint function
rome_constraint(obj);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
