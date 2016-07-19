function val = isempty(obj)

% PROF_VAR\ISEMPTY Overrides builtin isempty functionality by checking
% whether the AffineMap of obj is an empty matrix
%
%   is_empty_flag = isempty(A);
%
%
% Modification History: 
% 1. Joel 

val = isempty(obj.BiAffineMap);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
