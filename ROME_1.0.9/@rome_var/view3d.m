function view3d(obj)

% ROME_VAR\VIEW3D Development utility function which views the BiAffinemap
% of the current object in 3D
%
% Modification History: 
% 1. Joel 

% number of cols in each sub-block
inner_step = obj.NumMappedVars + 1;

% get row and cols of non-zeros
[r, c] = find (obj.BiAffineMap);

% z-coordinate
z = ceil(c ./ inner_step);

% change c
actual_c = mod(c - 1, inner_step) + 1;

% scatter
figure;
h = stem3(r, actual_c, z, 'filled');
xlabel('Out');
ylabel('Certain');
zlabel('Uncertain');


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
