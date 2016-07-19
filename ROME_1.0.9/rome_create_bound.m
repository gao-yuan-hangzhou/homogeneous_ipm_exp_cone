function bound_obj = rome_create_bound(y, varargin)

% ROME_CREATE_BOUND constructs a unified upper bound for E(y^+) from a list of bounds
% using infimal convolution. 
%
% USAGE:
% w = rome_create_bound(y, @rome_covar_bound);
% w = rome_create_bound(y, @rome_covar_bound, @rome_dirdev_bound);
% w = rome_create_bound(y, @rome_covar_bound, @rome_supp_bound);
% w = rome_create_bound(y, @rome_supp_bound, @rome_covar_bound, @rome_dirdev_bound);
%
% Input: y - LDR variable to be bounded
% Input: varargin - list of function handles to bounds
%
% Output: bound_obj - deterministic rome_var object constructed from bounds
% 
% Modification History: 
% 1. Joel 

h = rome_get_current_model();

% by default, if the varargin is empty, check the model to see which bounds
% to use
if(isempty(varargin) || isempty(varargin{1}))
    % intialize
    varargin = {};

    % some support info
    if(~(isempty(h.rndLB) && isempty(h.rndUB) ...
       && isempty(h.rndLC) && isempty(h.rndQC)))
        varargin = {@rome_supp_bound};
    end
    
    % some covar info
    if(~isempty(h.rndCovar))
        varargin = {varargin{:}, @rome_covar_bound};
    end
    
    % some dirdev info
    if(~(isempty(h.rndFDev) && isempty(h.rndBDev)))
        varargin = {varargin{:}, @rome_dirdev_bound};
    end
end

% initialize bounding object and sum of y's
bound_obj = 0; 
sum_y = 0;

z = h.get_rand_vars();

% iterate over the bounds to be combined
for ii = 1:numel(varargin)-1
     cur_bound = varargin{ii};
     aux_y = rome_linearrule(size(y), z);
     bound_obj = bound_obj + cur_bound(aux_y);
     sum_y = sum_y + aux_y;
end
 
% make the final bound
last_bound = varargin{end};
bound_obj = bound_obj + last_bound(y - sum_y);

% % TEMP:
% % iterate over the bounds to be combined
% for ii = 1:numel(varargin)
%      cur_bound = varargin{ii};
%      aux_y = rome_clone(y);         % make a new variable
%      bound_obj = bound_obj + cur_bound(aux_y);  % apply bound
%      sum_y = sum_y + aux_y;
% end
%  
% % make the constraint
% rome_constraint(sum_y == y);



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
