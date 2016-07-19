% ROME_CONSTANTS Describes the set of constants in ROME
%
%
% Modification History: 
% 1. Changed 'random' key string to 'uncertain' (Joel)
% 

classdef rome_constants
    % CONSTANTS
    % ----------
    properties (Constant = true)
        % Constants representing Conic Constraints
        NO_CONE  = -1;  % No Cone Restrictions
        ZERO     =  0;  % Restricted to 0
        NNOC     =  1;  % Nonnegative Cone
        SOC      =  2;  % Second Order Cone
        
        % Constant representing Continuity Constraints
        CONTINUOUS  = 0; % Continuous Variable
        INTEGER     = 1; % Integer Variable
        BINARY      = 2; % Binary Variable
        
        % list of string constants
        CONE_CONSTANTS = 1;
        ConeStr  = {'nonneg', 'lorentz'};
        
        CONTINUITY_CONSTANTS = 2;
        ContinuityStr  = {'integer', 'binary'};
        
        VARIABLE_CONSTANTS = 3;
        VariableStr  = {'uncertain', 'linearrule', 'empty'};
        
        % all constants
        AllStr = cat(2, rome_constants.ConeStr, ...
                        rome_constants.ContinuityStr, ...
                        rome_constants.VariableStr);
        
        EMPTY      = -1;
        CERTAIN    = 0;
        UNCERTAIN  = 1;
        LINEARRULE = 2;        
    end
    
    methods(Static = true)
        % checks if a given string represents a valid constant
        function valid_flag = is_valid_str(str)
            valid_flag  = any(strcmpi(str, rome_constants.AllStr));
        end
        
        % TRANSLATE: translates string into input format for rome_var
        % constructors
        function [var_type, str_options, arg_val] = translate(attrib_arr)
            % default: no args and certain variable type
            str_options = {};
            arg_val = {};
            set_cone_flag = false;
            set_continuity_flag = false;
            set_var_flag = false;
            
            % find variable type first
            var_type = rome_constants.CERTAIN;
            for ii = 1:numel(rome_constants.VariableStr)
                if(any(strcmpi(rome_constants.VariableStr{ii}, attrib_arr)))
                    if(set_var_flag)
                         error('rome_constants:translate:DuplicateAttributes', ...
                            'More than 1 Variable Type Attribute Specified');
                    end
                    
                    % special case: empty variable
                    if(ii == numel(rome_constants.VariableStr))
                        var_type = rome_constants.EMPTY;
                        return;
                    else                    
                        var_type = ii;                        
                    end
                    set_var_flag = true;
                end
            end
            
            % iterate over all the supplied attributes
            for ii = 1:numel(attrib_arr)
                % define str
                str = attrib_arr{ii};
               
                % find cone
                ind = find(strcmpi(str, rome_constants.ConeStr));
                if(~isempty(ind))
                    if(set_cone_flag)
                        error('rome_constants:translate:DuplicateAttributes', ...
                                'More than 1 Cone Attribute Specified');
                    end
                    if(set_var_flag && ind == rome_constants.SOC)
                        error('rome_constants:translate:NotYetImplement', ...
                            'Not Yet implement Lorentz constraint for uncertain and linearrrule variables');
                    end
                    % str_options = [str_options, ', ''Cone'', ', num2str(ind)];
                    str_options = cat(2, str_options, 'Cone'); 
                    arg_val = cat(2, arg_val, ind);
                    set_cone_flag = true;
                end
                
                % find continuity
                ind = find(strcmpi(str, rome_constants.ContinuityStr));
                if(~isempty(ind))
                    if(set_continuity_flag)
                        error('rome_constants:translate:DuplicateAttributes', ...
                            'More than 1 Continuity Attribute Specified');
                    end
                    if(set_var_flag)
                        error('rome_constants:translate:NoContinuity', ...
                            '"%s" attribute is not allowed for uncertain and linearrule variables', ...
                            rome_constants.ContinuityStr{ind});
                    end
                    % str_options = [str_options, ', ''Continuity'', ', num2str(ind)];
                    str_options = cat(2, str_options, 'Continuity');
                    arg_val = cat(2, arg_val, ind);
                    set_continuity_flag = true;
                end
            end
        end
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
