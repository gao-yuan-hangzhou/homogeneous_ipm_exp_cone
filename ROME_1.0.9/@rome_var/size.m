% ROME_VAR\SIZE Implements size method for rome_var objects
%
%
% Modification History: 
% 1. Joel 

function sz = size(obj, index)
    if(nargin == 1)
        sz = obj.Size;
    else
        sz = obj.Size(index);
    end    
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
