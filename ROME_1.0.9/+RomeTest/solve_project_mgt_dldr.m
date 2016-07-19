function [obj_val, x_val, y_val, r_val] = solve_project_mgt_dldr(C, beta, W, H, a, b, c)

% 
% +ROMETEST\SOLVE_PROJECT_MGT_DLDR Helper routine to solve
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

% number of graph edges
num_edges = (H - 1) * W + (W - 1) * H;

% begin rome environment
rome_begin;
h = rome_model('Project Management Example - DLDR');

% Define bounds of primitive uncertainties
z_LB = repmat(-1.2 / (2 * (1-beta)), num_edges, 1);
z_UB = repmat( 1.2 / (2 * beta)    , num_edges, 1);

% define covariance
z_cov = speye(num_edges) ./ (4 * beta * (1 - beta));
% z_std = speye(num_edges) ./ (2 * sqrt(beta * (1 - beta)));

% uncertain vars corresponding to edges in graph
z = rome_rand_model_var(num_edges);  
rome_constraint(z >= z_LB);
rome_constraint(z <= z_UB);
z.Covar = z_cov;
z.set_mean(0);


% model variables
x = rome_model_var(num_edges, 'Cone', rome_constants.NNOC);
rome_constraint(x <= 1);

y = rome_linearrule([H, W], z);
rome_constraint(y(1, 1) == 0);

% budget constraint
rome_constraint(sum(c.*x(:)) <= C );

% DLDR constraint
% g = rome_model_var(num_edges);
% [r, g] = RobustBounds.dldr_bound(num_edges, z, z_LB, z_UB);
% [r, g] = RobustBounds.dldr_bound(num_edges, z);

% test
r = rome_linearrule(num_edges, z);
% sel_ind = 1:3;
% I = num_edges;
% g = RobustBounds.meanpositivebound(-r, z, I, sel_ind); % currently not so efficient (but it works... so what's wrong?
g = rome_create_bound(-r, @rome_supp_bound, @rome_covar_bound);

% seqential constraint
dy1 = diff(y, [], 1); 
dy2 = diff(y, [], 2);
rome_constraint([dy1(:); dy2(:)] - r == b + a * (1 - x) .* z);

% objective
rome_minimize(strip_rand(y(H, W)) + sum(g));

% solve
h.solve;

% define output values
obj_val = h.objective;
x_val   = h.eval(x);
y_val   = linearpart(h.eval(y));
r_val   = linearpart(h.eval(r));

rome_end;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
