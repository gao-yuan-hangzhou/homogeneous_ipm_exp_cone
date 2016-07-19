 function obj = reshape(obj, varargin)

% ROME_VAR\RESHAPE reshapes the object to fit a new size
%
% reshaped_obj = reshape(obj, new_sz)
% reshaped_obj = reshape(obj, N1, N2, N3... )
%
% Modification History: 
% 1. Joel 

% parse arguments
if(numel(varargin) == 1)
    new_sz = varargin{1};
else
    new_sz = horzcat(varargin{:});
end

% error checking
if(prod(new_sz) ~= obj.TotalSize)
    error('rome_var:reshape:InvalidSize', 'Object must be reshaped to the same size');
end

if(any(int64(new_sz) ~= new_sz))
    error('rome_var:reshape:ReqIntegerSize', 'Size must be real integer values');
end

% squeeze out trailing singletons
while((length(new_sz) > 2) && (new_sz(end) == 1))
    new_sz = new_sz(1:end-1);
end

% if all ok, change the Size of the object
obj.Size = new_sz;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
