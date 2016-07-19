function bound_obj = rome_pwl_create_bound(y, b_vec, a_vec, varargin)

% ROME_PWL_CREATE_BOUND constructs a unified upper bound for piecewise
% linear (pwl) function of LDR
%
% Input: y - LDR variable to be bounded
% Input: varargin - list of function handles to bounds
%
% Output: bound_obj - deterministic rome_var object constructed from bounds
% 
% Modification History: 
% 1. Joel 

h = rome_get_current_model();

% check that y is column vector
if(~isvector(y) && size(y, 2) ~= 1)
    error('rome_pwl_create_bound:OnlyColVec', 'y must be a column vector');
end

% by default, if the varargin is empty, check the model to see which bounds
% to use
if(isempty(varargin) || isempty(varargin{1}))
    % intialize
    varargin = {};

    % some support info
    if(~(isempty(h.rndLB) && isempty(h.rndUB) ...
       && isempty(h.rndLC) && isempty(h.rndQC)))
        varargin = {@rome_pwl_supp_bound};
    end
    
    % some covar info
    if(~isempty(h.rndCovar))
        varargin = [varargin, {@rome_pwl_covar_bound}];
        % varargin = {varargin{:}, @rome_pwl_covar_bound};
    end
    
    % some dirdev info
    if(~(isempty(h.rndFDev) && isempty(h.rndBDev)))
        error('Not Implemented Yet');
        % varargin = {varargin{:}, @rome_dirdev_bound};
    end
end

% initialize bounding object and sum of y's
z = h.get_rand_vars();

% % iterate over the bounds to be combined
% for ii = 1:numel(varargin)-1
%      cur_bound = varargin{ii};
%      aux_y = rome_linearrule(size(y), z);
%      aux_b = rome_model_var(size(b_vec));
%      
%      bound_obj = bound_obj + cur_bound(aux_y, aux_b, a_vec);
%      sum_y = sum_y + aux_y;
%      sum_b = sum_b + aux_b;
% end
%  
% % make the final bound
% last_bound = varargin{end};
% bound_obj = bound_obj + last_bound(y - sum_y, b_vec - sum_b, a_vec);

% check output dimension
is_b_condensed = isvector(b_vec);

% check if a_vec is a vector, preprocess
% checks if a is specified in full or condensed format
is_a_condensed = isvector(a_vec); 
szA = size(a_vec);
if(~is_a_condensed)    
    % reshaped_a_vec for use when a is fully specified
    reshaped_a_vec = reshape(a_vec, [szA(1), szA(2) * szA(3)]);
    reshaped_a_vec = blk_transpose(reshaped_a_vec, szA(1), szA(2));    
end

if(~is_b_condensed)
    out_dim = size(b_vec, 1);
elseif(~is_a_condensed)
    out_dim = size(a_vec, 1);
else
    out_dim = size(y, 1);
end


% iterate over the bounds to be combined
bound_obj = 0; 
sum_C = 0;
sum_Y = 0;

% problem with a_vec!! Stopped here!
for ii = 1:numel(varargin)
    cur_bound = varargin{ii};
    aux_y = rome_linearrule(size(y), z);         % make a new variable
    aux_b = rome_model_var(size(b_vec));
    if(is_b_condensed)
        aux_b = repmat(aux_b, [out_dim, 1]);
    end
   
    bound_obj = bound_obj + cur_bound(aux_y, aux_b, a_vec);  % apply bound

    sum_Y = sum_Y + aux_y;
    sum_C = sum_C + aux_b;
%     sum_Y = sum_Y + strip_certain(aux_y);
        
%     if(is_a_condensed)
%         % optimization for case when entries of a are scalars
%         sum_C = sum_C + strip_rand(aux_y) * a_vec + aux_b;
%     else       
%         sum_C = sum_C + reshape(reshaped_a_vec * strip_rand(aux_y), out_dim, szA(3)) + aux_b;
%     end
end
 
% preprocess before constraints
if(is_b_condensed)
    b_vec = repmat(b_vec, [out_dim, 1]);
end

% if(is_a_condensed)
%     rome_constraint(sum_C == strip_rand(y) * a_vec + b_vec);
% else
%     rome_constraint(sum_C == reshape(reshaped_a_vec * strip_rand(y), out_dim, szA(3)) + b_vec);
% end

% % make the constraints
% rome_constraint(sum_Y == strip_certain(y));

rome_constraint(sum_Y == y);
rome_constraint(sum_C == b_vec);

    
% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.