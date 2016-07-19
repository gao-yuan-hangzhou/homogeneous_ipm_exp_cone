function bDef = isdeflected(sol_obj)
%
% ROME_SOL\ISDEFLECTED Returns true if the current object has deflections,
% otherwise returns false
%
% Usage:
% bIsDeflected = isdeflected(sol_obj);
%
% History
% 1. Created by Joel 19 May 2009 
%

bDef = ~isempty(find(sol_obj.DeflectCoeff, 1));


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
