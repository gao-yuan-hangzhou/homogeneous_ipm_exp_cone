function obj = strip_certain(obj)

% PROF_VAR\STRIP_CERTAIN Retains only the uncertain part of object, though as
% a deterministic variable. 
%
% Note that this returns the transpose: i.e. for a vector-valued LDR,
%
% y(z) = y_0 + Y * z
%
% This returns Y' instead.
%
%
% Modification History: 
% 1. Joel 

blk_size =(obj.NumMappedVars+1);
obj.BiAffineMap = blk_transpose(obj.BiAffineMap, obj.TotalSize, blk_size); % outer transpose
obj.BiAffineMap = obj.BiAffineMap((obj.TotalSize+1):end, :);

% modify object's size
obj.Size = [obj.TotalSize, obj.NumMappedRandVars]; 
obj.DiagMult = NaN;
obj.NumUnmappedRandVars = 0;
obj.NumMappedRandVars = 0;

% return the transpose of the object
obj = obj';


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
