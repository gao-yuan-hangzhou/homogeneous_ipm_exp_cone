function [deflect_coeff, deflect_val] = deflectedpart(sol_obj)
%
% ROME_SOL\DEFLECTEDPART Returns the deflected part of the solution object.
% Returns both the coefficients and the deflected value
%
% Usage:
% [x_def_coeff, x_def_vals] = deflectedpart(x);
%
% History
% 1. Created by Joel 19 May 2009 
%

deflect_coeff = sol_obj.DeflectCoeff;
deflect_val = sol_obj.DeflectAffineMap;

% reshape deflect_coeff to size
new_sz = [sol_obj.Size, size(deflect_coeff, 2)];
deflect_coeff = reshape(deflect_coeff, new_sz);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
