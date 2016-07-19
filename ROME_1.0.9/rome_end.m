function rome_end()

% ROME_END Clears the rome environment for good. To be called to free up
% any existing rome resources

global ROME_ENV;
if(~isempty(ROME_ENV))    
    for ii = 1:length(ROME_ENV.arr_models)
        delete(ROME_ENV.arr_models(ii));
    end

    ROME_ENV.arr_models = [];
    ROME_ENV.curr_model = [];
    ROME_ENV.num_models = 0;
end
% clear -GLOBAL ROME_ENV;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
