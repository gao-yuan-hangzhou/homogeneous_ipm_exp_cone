function varargout = rome_begin(varargin)

% ROME_BEGIN Starts up the ROME Environment. Needs to be called before any
% other function in ROME.

global ROME_ENV;

% clear any existing models
if(~isempty(ROME_ENV))
    if(isfield(ROME_ENV, 'arr_models'))
        for ii= 1:length(ROME_ENV.arr_models) 
            delete(ROME_ENV.arr_models(ii));
        end
    end
end

ROME_ENV.arr_models = [];
ROME_ENV.curr_model = [];
ROME_ENV.num_models = 0;

% only create on first attempt
if(~isfield(ROME_ENV, 'VERBOSE'))
    ROME_ENV.VERBOSE    = 1;
end

% make a structure of parameters to supply to underlying solvers
% only create on first attempt
if(~isfield(ROME_ENV, 'DEF_SOLVER'))
    ROME_ENV.DEF_SOLVER = 'CPLEX';    % Default Solver is CPLEX Solver
end

% make a structure of parameters to supply to underlying solvers
% only create on first attempt
if(~isfield(ROME_ENV, 'MSK_PARAMS'))
    ROME_ENV.MSK_PARAMS = []; % MOSEK Params;
end

% if argument is supplied
if(nargout == 0)
    return;
elseif(nargout == 1)
    varargout{1} = rome_model(varargin{:});
else
    error('rome_begin:WrongNumOutArgs', 'Too many output arguments');    
end

% notice that this does not actually delete the memory, but
% simple hides it from the current workspace
clear ROME_ENV;

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
