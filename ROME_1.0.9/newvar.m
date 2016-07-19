function out_obj = newvar(varargin)

%
% NEWVAR Global Function to create model-level variables in rome. All the 
% variables are created in the calling function's workspace.
% 
% A) Vanilla: 
%   1. NEWVAR('x') creates a scalar variable x
%   2. NEWVAR('x(N)') creates a (column) vector variable x of size N
%   3. NEWVAR('x(N1, N2, .. )') creates an N-D variable of size N1 x N2 ..
% 
% B) With Attributes: 
%   1. NEWVAR('x', 'nonneg') creates a non-negative variable x
%   2. NEWVAR('x', 'binary') creates a binary variable x
%   3. NEWVAR('x', 'integer') creates an integer variable x
%
% Attributes can be combined together and used with size information to
% form more complex variable creations: e.g.
%   4. NEWVAR('x(3, 3)', 'nonneg', 'integer') creates a 3x3 matrix of
%      integer variables, which are constrained to be non-negative.
% 
% See the User's Guide for a complete list of allowable attributes.
%
% C) Uncertain Variables:
%   1. NEWVAR('z', 'uncertain') creates a scalar uncertain variable z.
%   2. NEWVAR('z(N)', 'uncertain') creates a vector uncertain variable z, size N
%   3. NEWVAR('z(N1, N2)', 'uncertain') creates a N-D uncertain variable z of size N1 x N2 ...
% 
% Uncertain variables can have conic restriction attributes, but not
% continuity attributes.
%   4. NEWVAR('x(3)', 'uncertain', 'nonneg') creates a 3x1 vector of
%      uncertain variables, which are constrained to be non-negative.
%   5. NEWVAR('x(10)', 'uncertain', 'lorentz') creates a 10x1 vector of
%      uncertain variables, which is confined to the Lorentz-cone (Second
%      Order Cone). (NOT IMPLEMENTED YET)
% 
% D) Linear Decision Rules
%   1. NEWVAR('x(z)', 'linearrule') creates a scalar linearrule variable 
%      x, which depends on z.
%   2. NEWVAR('x(N, z)', 'linearrule') creates a length N vector linearrule
%      variable x, which depends on z.
%   3. NEWVAR('x([N1, N2.. ], z)', 'linearrule') creates a N-D linearrule 
%      variable z of size N1 x N2 ...
%   4. NEWVAR('x(N1, N2.. , z)', 'linearrule') also creates a N-D linearrule 
%      variable z of size N1 x N2 ...
% 
% Linearrule variables can have conic restriction attributes but not
% continuity attributes. 
%   4. NEWVAR('x(3, z)', 'linearrule', 'nonneg') creates a 3x1 vector of
%      uncertain variables, which are constrained to be non-negative.
%   5. NEWVAR('x(10, z)', 'linearrule', 'lorentz') creates a 10x1 vector of
%      uncertain variables, which is confined to the Lorentz-cone (Second
%      Order Cone). (NOT IMPLEMENTED YET)
% 
% Notice for linearrules, the uncertain variable z MUST be defined in the
% calling function's workspace BEFORE x(z) can be called.
% 
% E) Empty Variables: 
%   1. NEWVAR('x(N1, N2, .. )', 'empty') creates an empty N-D variable,
%   which can be used as a container. If other attributes are present, they
%   will be ignored.
%
% General Comments: Attribute strings are case insensitive. Thus 'UNCERTAIN',
% 'uncertain', and 'UnCeRTaiN', are all equally acceptable strings;
%
% Modified by:
% 1. Joel 

% check which ones are valid strings
valid_arr = cellfun(@rome_constants.is_valid_str, varargin);

% default vals
attrib_arr = {};
new_var_arr = varargin;

% split array into attributes and variables
first_attrib = find(valid_arr, 1);
if(~isempty(first_attrib))
    attrib_arr = varargin(first_attrib:end);
    new_var_arr = varargin(1:first_attrib-1);
end

% check that constants are in 1 contiguous block
if(~all(valid_arr(first_attrib:end)))
    errant_index = first_attrib-1 + find(~valid_arr(first_attrib:end), 1);
    error('new_var:InvalidConstants', '"%s" not a valid attribute string', ...
            varargin{errant_index});
end

% if all ok, make arguments into the rome_model_var function
[var_type, str_options, arg_vals] = rome_constants.translate(attrib_arr);

% switch between creation functions
switch(var_type)
    case rome_constants.CERTAIN
        fn_name = 'rome_model_var';
    case rome_constants.UNCERTAIN
        fn_name = 'rome_rand_model_var';
    case rome_constants.LINEARRULE
        fn_name = 'rome_linearrule';
    case rome_constants.EMPTY
        fn_name = 'rome_empty_var';
    otherwise
        error('new_var:UnrecognizedType', 'Unrecognized creation type');
end

% find numeric arguments (if there are any, we are in function mode)
if(isempty(new_var_arr))
    new_var_arr = {1};
end
is_numeric_arr = cellfun(@isnumeric, new_var_arr);
is_rome_var_arr= cellfun(@(x) isa(x, 'rome_var'), new_var_arr);
if(any(is_numeric_arr) || any(is_rome_var_arr))
    if(~all(is_numeric_arr | is_rome_var_arr))
        error('new_var:IncorrectInput', 'Need numeric size arguments in new_var function mode');
    end
    
    % check for rome var and find index, replace actual var with name
    opt_array = cat(1, str_options, arg_vals);
    sz_args = [new_var_arr{~is_rome_var_arr}];
    if(any(is_rome_var_arr))
        dep_rand_var = new_var_arr{is_rome_var_arr};
        out_obj = feval(fn_name, sz_args, dep_rand_var, opt_array{:});
    else
        out_obj = feval(fn_name, sz_args, opt_array{:});
    end    
else
    % make options
    str_options = arrayfun(@(x) [', ''', x{:}, ''''], str_options, 'UniformOutput', false);
    arg_vals    = arrayfun(@(x) [', ', num2str(x{:})], arg_vals, 'UniformOutput', false);
    str_options = cat(1, str_options, arg_vals);
    str_options = [str_options{:}];
    
    % iterate over new arguments
    for ii = 1:numel(new_var_arr)
        % find parentheses
        str_size = '1, 1';
        cur_var = new_var_arr{ii};
        start_ind = min(findstr(cur_var, '('));
        end_ind = max(findstr(cur_var, ')'));

        % check the start and end indices of name string
        if(~isempty(start_ind) && ~isempty(end_ind))
            str_size = cur_var(start_ind+1:end_ind-1);
            cur_var  = cur_var(1:start_ind-1);
        end
        
        % make command string
        cmd_str = [fn_name, '(', str_size, str_options, ')'];

        % evalate
        assignin('caller', cur_var, evalin('caller', cmd_str));
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
