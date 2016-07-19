function val = mean(obj)

% PROF_VAR\MEAN Returns the mean of the uncertain variable
%
%
% Modified By:
% 1. Joel

% val = obj.Mean;

% Error check uncertainness 
if(~(obj.IsRand || obj.IsLDR))
    error('rome_var:mean:ReqRandLDR', 'Object should be uncertain or LDR to get mean');
end

% get from global variable
h = rome_get_current_model();
primitive_mean = h.get_mean(obj);

% check if mean is constant
[r, c] = find(primitive_mean.BiAffineMap);
if(all(c == 1))
    % if all col indices are 1, then it's a constant
    primitive_mean = full(primitive_mean.BiAffineMap(:, 1));
end

if(obj.IsRand)
    % insert value and return
    val = obj.insert(primitive_mean);
else
    val = strip_rand(obj) + reshape(strip_certain(obj)' * primitive_mean, size(obj));
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
