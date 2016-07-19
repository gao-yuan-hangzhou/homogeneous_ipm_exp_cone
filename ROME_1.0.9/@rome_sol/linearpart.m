function val = linearpart(sol_obj)
%
% ROME_SOL\LINEARPART Returns the linearpart of the solution object
%
% Usage:
% x_val = x.linearpart();
% x_val = linearpart(x);
%
% History
% 1. Created by Joel 19 May 2009 
%

val = sol_obj.LDRAffineMap;

% reshape object to size
new_sz = [sol_obj.Size, size(val, 2)];
val = reshape(val, new_sz);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
