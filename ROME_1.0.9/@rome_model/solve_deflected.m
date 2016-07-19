function [x_min, f_min, sol_stat, details] = solve_deflected(model_obj, varargin)
%
% ROME_MODEL\SOLVE_DEFLECTED Performs Bi-deflection and solve in a single
% step
%
% Modified by:
% 1. Joel (Created 25 May 2009)

% check for maximization
bMaxFlag = false;
if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
    model_obj.ObjFn = -model_obj.ObjFn;
    model_obj.MinMaxFlag = rome_model.MINIMIZE;
    bMaxFlag = true;
end

% use deflection options from model object. 
% apply_na_bdldr has in-built handlers for default options
model_obj.apply_na_bdldr(model_obj.defObjFn, ...        % subproblem objective 
                         model_obj.defLDRRestrict, ...  % LDR Restriction
                         model_obj.defBounds);          % Array of bounds

% call solve
[x_min, f_min, sol_stat, details] = model_obj.solve(varargin{:});

if(bMaxFlag)
    f_min = -f_min;
    model_obj.ObjVal = -model_obj.ObjVal;
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
