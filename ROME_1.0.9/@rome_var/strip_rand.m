function obj = strip_rand(obj)

% PROF_VAR\STRIP_RAND Retains only the deterministic part of object
%
% TEMP!
%
% Modification History: 
% 1. Joel 

obj.BiAffineMap = obj.BiAffineMap(:, 1:(obj.NumMappedVars+1));
obj.NumUnmappedRandVars = 0;
obj.NumMappedRandVars = 0;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
