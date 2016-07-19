function sz = numdepvars(sol_obj)
%
% ROME_SOL\NUMDEPVARS Returns Number of Dependent Variables
%
% Usage: 
% sz = size(x);
%
% History
% 1. Created by Joel 14 May 2009 
%

sz = size(sol_obj.LDRAffineMap, 2) - 1;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
