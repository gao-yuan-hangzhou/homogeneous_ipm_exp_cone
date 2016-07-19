function h_model_object = rome_get_current_model()

% ROME_GET_CURRENT_MODEL Returns a handle to the currenly selected model
% object.
%
% Will raise an exception if the ROME Environment has not been initialized
% with rome_begin.
%
% Modification History: 
% 1. Joel 

global ROME_ENV;

try
    % try to assign model object    
    h_model_object = ROME_ENV.curr_model;
catch ME
    
    % if environment hasn't been initialized, give a more descriptive error
    % and exit
    if(strcmp(ME.identifier, 'MATLAB:nonStrucReference'))
        error('rome_get_current_model:MissingGlobal', ...
            'Global ROME Environment not started yet. Have you called rome_begin?');
    else
        rethrow(ME);
    end
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
