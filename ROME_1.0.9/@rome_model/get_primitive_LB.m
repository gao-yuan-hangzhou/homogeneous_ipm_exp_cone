function bnd = get_primitive_LB(model_obj, var_obj)

% PROF_MODEL\GET_PRIMITIVE_LB Returns the lower bound of the PRIMITIVE
% uncertainty for the supplied variable object.
%
% 
% Modified by: 
% 1. Joel Goh
%

start_ind = var_obj.NumUnmappedRandVars + 1;
end_ind   = var_obj.NumUnmappedRandVars + var_obj.NumMappedRandVars;

% expand vector if necessary
num_append = pos(end_ind - numel(model_obj.rndLB));
if(num_append > 0)
    model_obj.rndLB = [model_obj.rndLB; -Inf(num_append, 1)];
end

% return bound here
bnd = model_obj.rndLB(start_ind:end_ind);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
