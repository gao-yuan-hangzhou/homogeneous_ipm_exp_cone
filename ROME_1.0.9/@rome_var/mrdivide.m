function out_var_obj = mrdivide (A, B)

% ROME_VAR\MRDIVIDE Implements matrix (right) standard division in ROME
%
%   C = A / B 
%
%   A must be a rome_var
%   B must be numeric
%   A, B must be of the same size. 
%   Currently mirrors A ./ B
%
% Modification History: 
% 1. Joel 

out_var_obj = A ./ B;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
