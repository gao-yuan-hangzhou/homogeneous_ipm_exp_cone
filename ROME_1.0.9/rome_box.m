function rome_box(obj, lb, ub)

% ROME_BOX Applies a box-type constraint, lb <= obj <= ub on the variable
%
%   obj must be a scalar-valued rome_var or a real scalar. 
%
% Modification History: 
% 1. Joel 

rome_constraint(obj >= lb);
rome_constraint(obj <= ub);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
