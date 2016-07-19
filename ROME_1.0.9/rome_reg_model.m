function rome_reg_model(h_model_object)

% ROME_REG_MODEL
% REGISTERS THE CURRENT MODEL IN THE GLOBAL ENVIRONMENT


global ROME_ENV;

ROME_ENV.arr_models = [ROME_ENV.arr_models, h_model_object];
ROME_ENV.curr_model = h_model_object;
ROME_ENV.num_models = ROME_ENV.num_models + 1;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
