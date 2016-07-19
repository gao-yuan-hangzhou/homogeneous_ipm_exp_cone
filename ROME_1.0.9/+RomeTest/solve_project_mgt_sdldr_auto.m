function [obj_val, x_val, y_val] = solve_project_mgt_sdldr_auto(C, beta, W, H, a, b, c, z_pos_mean)

% 
% +ROMETEST\SOLVE_PROJECT_MGT_SDLDR_AUTO Helper routine to solve
% project management example using a segregated-deflected linear decision rule
%
%   [obj_val, x_val, y_val] = SOLVE_PROJECT_MGT_SDLDR_AUTO(C, beta, W, H, a, b, c)
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
%   z_pos_mean   : Mean of positive part of primitive uncertainty
%
%   returned x_val will be a W by H matrix, and y_val will be a NUM_EDGES
%   by NUM_EDGES matrix, where NUM_EDGES = (H-1)*W + (W-1)*H.
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

% define default parameters
if(nargin < 8)
    z_pos_mean = 0.5; % mean of positive part of primitive uncertainty
end
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
h = rome_model('Project Management Example - SDLDR (AUTO)');

% define number of edges 
num_edges = (H - 1) * W + (W - 1) * H;

% Define bounds of primitive uncertainties
z_LB = -1.2 / (2 * (1-beta));
z_UB =  1.2 / (2 * beta);
z_seg_LB = [repmat(-z_pos_mean, num_edges, 1); repmat(z_LB + z_pos_mean, num_edges, 1)];
z_seg_UB = [repmat(z_UB - z_pos_mean, num_edges, 1); repmat(z_pos_mean, num_edges, 1)];

% Compute Covariance Matrix 
z_var_1  = 0.25 * ((1 - beta) ./ beta);
z_var_2  = 0.25 * (beta ./ (1 - beta));
z_cross  = z_pos_mean.^2;
z_seg_cov = [speye(num_edges) * z_var_1, speye(num_edges) * z_cross; ...
             speye(num_edges) * z_cross , speye(num_edges) * z_var_2];
z_seg_cov = z_seg_cov + speye(2*num_edges)*eps;
% z_seg_std    = chol(z_seg_cov); % Factorize Covariance Matrix

% uncertain vars corresponding to edges in graph
z = rome_rand_model_var(2*num_edges);
z1 = z(1:num_edges);
z2 = z(num_edges+1:end);
rome_constraint(z >= z_seg_LB);
rome_constraint(z <= z_seg_UB);
z.Covar = z_seg_cov;
z.set_mean(0);

% model variables
x = rome_model_var(num_edges, 'Cone', rome_constants.NNOC);
rome_constraint(x <= 1);

y = rome_linearrule([H, W], z);
rome_constraint(y(1, 1) == 0);

% seqential constraint
dy1 = diff(y, [], 1); 
dy2 = diff(y, [], 2);
rome_constraint([dy1(:); dy2(:)] >= b + a * (1 - x) .* (z1 + z2));

% budget constraint
rome_constraint(sum(c.*x(:)) <= C );

% objective
rome_minimize(strip_rand(y(H, W)));  %+ sum(strip_rand(w)));

% make dldr
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
