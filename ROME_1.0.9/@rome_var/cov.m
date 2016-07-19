function [proj_covar, covar_mix] = cov(obj)

% PROF_VAR\COV Returns the covariance of the uncertain variable
%
%
% Modified By:
% 1. Joel

h = rome_get_current_model();
[proj_covar, covar_mix] = h.get_cov(obj);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
