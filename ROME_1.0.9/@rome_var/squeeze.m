function obj = squeeze(obj)

% PROF_VAR\SQUEEZE Removes Singleton dimensions
%
% TEMP
%
% Modification History: 
% 1. Joel 
% 2. Melvyn (18/09/2009) Fixed bug obj.Size is a scalar

if(length(obj.Size) <= 2)
    return;
else
    obj.Size(obj.Size == 1) = [];
    % Corrected by Melvyn
    if length(obj.Size) == 1;
        obj.Size = [obj.Size 1];
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
