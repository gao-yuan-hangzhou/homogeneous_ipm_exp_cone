function val = insert(sol_obj, z_val, vec_flag)
%
% ROME_SOL\INSERT Instantiates the current solution with uncertainties.
% Accepts a final argument, which controls whether we assume the input is
% vectorized. 
%
% USAGE:
% Suppose x_sol depends the 5 x 1 primitive uncertainties z, and has output
% dimension 3 x 1, and we want to instantiate it with 100 standard
% uniforms.
%
% % Case 1: No vectorize
% % --------------------
% x_val = zeros(3, 100);
% for ii = 1:100
%   x_val(:, ii) = x_sol.insert(rand(5, 1));
% end
%
% % Case 2: With vectorize
% % --------------c--------
% vectorize_flag = true;
% x_val = x_sol.insert(rand(5, 100), vectorize_flag);
%
% History
% 1. Created by Joel 14 May 2009 
%

if(nargin < 3)
    vec_flag = false;
end

if(~vec_flag)
    % Size check if we don't want to vectorize
    if(isvector(z_val))
        if(numel(z_val) ~= size(sol_obj.LDRAffineMap, 2) - 1)          % Mar 22 Fix by Joel
            error('rome_sol:insert:SizeMismatch', ...
                'Input uncertainties must have the same number of elements as rome_sol  object.');
        end
    else       
        if(ndims(z_val) ~= ndims(sol_obj.Size))
            error('rome_sol:insert:DimsMismatch', ...
                    'Input uncertainties must have same dimension as rome_sol object.');   
        end
        sz = size(z_val);
        if(prod(sz(1:end-1)) ~= size(sol_obj.LDRAffineMap, 2) - 1)      % Mar 22 Fix by Joel
            error('rome_sol:insert:SizeMismatch', ...
                    'Input uncertainties must have same size as rome_sol object.');
        end
    end
end

% vectorize input
OutDim = size(sol_obj.LDRAffineMap, 2) - 1;
N = ceil(numel(z_val) / OutDim);
z_val = [ones(1, N) ; reshape(z_val, [OutDim, N])];

% apply to output
if(sol_obj.isdeflected)
    val = sol_obj.LDRAffineMap * z_val + ... 
            sol_obj.DeflectCoeff * neg(sol_obj.DeflectAffineMap * z_val);
else
    val = sol_obj.LDRAffineMap * z_val;
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
