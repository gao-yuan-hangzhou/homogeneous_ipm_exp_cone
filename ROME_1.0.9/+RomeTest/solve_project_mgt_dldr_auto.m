function [obj_val, x_val, y_val] = solve_project_mgt_dldr_auto(C, beta, W, H, a, b, c)

% 
% +ROMETEST\SOLVE_PROJECT_MGT_DLDR_AUTO Helper routine to solve
% project management example using a linear decision rule
%
%   [obj_val, x_val, y_val] = SOLVE_PROJECT_MGT_LDR(C, beta, W, H, a, b, c)
%   returns the minimized expected project time, resources (x_val)
%   committed each node, and start times (y_val) of each node.
%
%   C   : Postive scalar value defining project budget
%   beta: Parameter used to determine bounds on primitive uncertainty
%   W   : Width of project grid (positive integer)
%   H   : Height of project grid (postive integer)
%   a   : Linear part of affine dependence of time on resources (scalar)
%   b   : Constant part of affine dependence of time on resources (scalar)
%   c   : Cost per unit of resources.
%
%   returned x_val will be a W by H matrix, and y_val will be a NUM_EDGES
%   by NUM_EDGES matrix, where NUM_EDGES = (H-1)*W + (W-1)*H.
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

% define default parameters
if(nargin < 7)
    c = 1;  % cost of each unit of resource (constant)
end
if(nargin < 6)
    b = 3;  % constant in relating time with resources
end
if(nargin < 5)
    a = 3;  % linear part for relating time with resources
end
if(nargin < 4)
    H = 4;  % Height of project grid
end
if(nargin < 3)
    W = 6;  % Width of project grid
end
if(nargin < 2)
    beta = 0.1; % parameter used to define uncertainty support
end
if(nargin < 1)
    C = 8;      % project budget 
end


% begin rome environment
rome_begin;
h = rome_model('Project Management Example - DLDR AUTO');

% Define bounds of primitive uncertainties
z_LB = -1.2 / (2 * (1-beta));
z_UB =  1.2 / (2 * beta);


% uncertain vars corresponding to edges in graph
num_edges = (H - 1) * W + (W - 1) * H;
z = rome_rand_model_var(num_edges);  
rome_box(z, z_LB, z_UB);

% make covariance matrix
z_cov = speye(num_edges) ./ (4 * beta * (1 - beta));

% assume covar
z.Covar = z_cov;
z.set_mean(0);


% model variables
x = rome_model_var(num_edges, 'Cone', rome_constants.NNOC);
rome_constraint(x <= 1);

y = rome_linearrule([H, W], z);
rome_constraint(y(1, 1) == 0);


% seqential constraint
dy1 = diff(y, [], 1); 
dy2 = diff(y, [], 2);
rome_constraint([dy1(:); dy2(:)] >= b + a * (1 - x) .* z);

% budget constraint
rome_constraint(sum(c.*x(:)) <= C );

% objective
rome_minimize(strip_rand(y(H, W)));

% call dldr here
% h.apply_dldr(strip_rand(y(H, W)), y, @rome_supp_bound, @rome_covar_bound);
h.apply_na_bdldr(strip_rand(y(H, W)), y, @rome_supp_bound, @rome_covar_bound);

% solve
h.solve;

% define output values
obj_val = h.objective;
x_val   = h.eval(x);
y_val   = linearpart(h.eval(y));

rome_end;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
