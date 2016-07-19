% Constructor
% ------------
function obj = rome_rand_model_var(varargin)

% ROME_MODEL_RAND_VAR Describes a model level variable in ROME
%
%
% Modification History:
% 1. Joel
%

% check that none of the arguments are rome_var objects. 
if(any(cellfun(@(x) isa(x, 'rome_var'), varargin)))
    error('rome_rand_model_var:NoVarInput', ...
        'rome_model_var does not accept rome_var inputs. Did you mean to call rome_linearrule instead?');
end

% construct the superclass
obj = rome_var(varargin{:});

% Get the current model and query its state
h_curr_model = rome_get_current_model();

try
    % define number of unmapped variables
%     obj.NumUnmappedVars = h_curr_model.NumVars;
    obj.NumUnmappedVars = 0;
    obj.NumUnmappedRandVars = h_curr_model.NumRandVars;
catch ME
    % Will Enter here if the user forgets to call rome_model
    % beforehand
    if(strcmp(ME.identifier, 'MATLAB:nonStrucReference'))
        warning('rome_model_var:ctor:UninitializedModel', ...
            'No Rome Model found: Creating a new model. To remove this error message, call rome_model first');
        h_curr_model = rome_model;
    else
        rethrow(ME);
    end
end

% make Uncertain variables
obj.NumMappedRandVars = obj.TotalSize;

% create the mapping rules
obj.BiAffineMap = [spalloc(obj.TotalSize, 1, 1), speye(obj.TotalSize)];

% All model variables are have identity mapping by default
obj.DiagMult = 1;

% add this variable to the current model
h_curr_model.add_new_var(obj);

% Allocate space for mean of uncertain variable
h_curr_model.ZRnd = [h_curr_model.ZRnd; zeros(obj.TotalSize, 1)];
h_curr_model.rndVarType = [h_curr_model.rndVarType; repmat('C', obj.TotalSize, 1)];


% mean is basically identical to original obj but with shifted indices
mean_obj = obj(:); 
mean_obj.NumUnmappedRandVars = obj.NumUnmappedRandVars + obj.TotalSize;

% % possibly extend mean
% h_curr_model.rndMean = [h_curr_model.rndMean; 
%                         zeros(pos(h_curr_model.NumRandVars - 2*obj.TotalSize - length(h_curr_model.rndMean)), 1)];
% 
% % construct pointer to mean
% if(isempty(h_curr_model.rndMean))
%     h_curr_model.rndMean = mean_obj;
% else
%     h_curr_model.rndMean = [h_curr_model.rndMean; mean_obj];
% end

% TESTING
h_curr_model.rndMean = [h_curr_model.rndMean; ...
                        mean_obj; zeros(obj.TotalSize, 1)];
% END TESTING

% attach mean flag
h_curr_model.ZIsMean = logical([h_curr_model.ZIsMean; false(obj.TotalSize,1); true(obj.TotalSize, 1)]);

% Apply the conic constraint at this time. Check for
% presence/absence of constraint done internally in
% rome_constraint function
rome_constraint(obj);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
