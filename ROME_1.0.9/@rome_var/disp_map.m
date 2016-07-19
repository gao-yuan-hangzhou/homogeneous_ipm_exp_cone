function A = disp_map(obj)
%ROME_VAR\DISP_MAP Display the biaffine map of a scalar rome_var object of
%size 1 (for debugging)

if (~isscalar(obj))
    error('rome_var:disp_map:NotScalar', ...
        'Object is not a scalar.');
end

if (obj.TotalSize ~= 1)
    error('rome_var:disp_map:TotalSizeIsNotOne', ...
        'Total size is not one.');
end

A = full(reshape(obj.BiAffineMap, ...
    obj.NumMappedVars+1, obj.NumMappedRandVars+1));
A = [ NaN, 0, obj.NumUnmappedRandVars+(1:obj.NumMappedRandVars); ...
    [0, obj.NumUnmappedVars+(1:obj.NumMappedVars)]', A];

end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
