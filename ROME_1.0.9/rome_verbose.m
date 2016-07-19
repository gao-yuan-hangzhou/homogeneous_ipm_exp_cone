function val = rome_verbose(varargin)

% ROME_VERBOSE Sets/Gets current verbosity level
% This controls the extend of output. 
%
% USAGE:
% 1. Setting verbosity level (returns old lvl)
%       lvl = 1
%       old_lvl = rome_verbose(lvl);
% 
% 2. Getting verbosity level:
%       lvl = rome_verbose();
%
% Current allowed levels (9 Mar 2009)
%   0 - Silent
%   1 - Outputs result of optimization
%
% Last Modified:
% 1. Joel 9 Mar 2009
% 

global ROME_ENV;

% check for empty
if(isempty(ROME_ENV))
    rome_begin;
end

% set return level
val = ROME_ENV.VERBOSE;

if(nargin > 1)
    error('rome_verbose:TooManyArguments', 'Cannot have more than 1 argument');
elseif(nargin > 0)
    ROME_ENV.VERBOSE = varargin{:};
end



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
