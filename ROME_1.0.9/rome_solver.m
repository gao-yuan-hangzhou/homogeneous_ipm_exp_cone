function varargout = rome_solver(varargin)

% ROME_SOLVER Sets / Gets current default solver
%
% USAGE:
% 1. Choosing a new solver:
%       new_solver = 'SDPT3DUAL'; % choose SDPT3 Dual solver
%       old_solver = rome_solver(new_solver); 
% 
% 2. Poll current solver:
%       old_solver = rome_solver();
%
% 3. Display a list of supported solver strings in a cell array
%       rome_solver('all');
%
%
% Last Modified:
% 1. Joel 9 Mar 2009
% 

% check if polling for all solvers
if(nargin == 1 && strcmp(varargin{1}, 'all'))
    fprintf('\nCurrently supported solver engines in ROME:\n');
    fprintf('\n--- Main ---\n');
    fprintf('CPLEXINT \nMOSEK \nSDPT3DUAL\n');
    fprintf('\n--- Experimental ---');
    fprintf('\nCPLEX \nSDPT3 \nCPLEXDUAL\n\n');    
    return;
end

% otherwise,
global ROME_ENV;

% check for empty
if(isempty(ROME_ENV))
    rome_begin;
end

% set return level
varargout{1} = ROME_ENV.DEF_SOLVER;
if(nargin > 1)
    error('rome_verbose:TooManyArguments', 'Cannot have more than 1 argument');
elseif(nargin > 0)
    ROME_ENV.DEF_SOLVER = varargin{1};    
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
