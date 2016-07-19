function len = length(obj)

% ROME_VAR\LENGTH Implements length method for rome_var objects
%
%
% Modification History: 
% 1. Joel 

if(isvector(obj))
    len = (obj.Size(1) * obj.Size(2));
else
    len = obj.Size(1);
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
