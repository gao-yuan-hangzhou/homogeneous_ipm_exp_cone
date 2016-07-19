% +ROMETEST\SIMPLE_WAREHOUSE Solves a robust warehouse problem using 
% a linear-decision rule.
%
% Reference: 
%
% Modified by:
% 1. Joel (22 Oct 2008)

% begin rome environment
rome_begin;

% default arguments
r_j  = [  1, 1.5,  2,  3, 100]';      % column-vector of retrieve costs
s_j  = [  1, 1.5,  2,  3, 100]';      % column-vector of store costs
c_j  = [ 15,  30, 45, 60, 300000]';   % column-vector of storage capacities
a_it = [35 35 35; ...   % fixed arrival of product i at time t
        35 35 35; ...
        10 10 25];  
mean_d_it = [29, 19, 45; ... % mean demand of product i at time t
             10, 47, 37; ...
             6,  7, 20];
z_UB = 1;
z_LB = -1;

% DATA
M = size(a_it, 1);  % number of products
N = length(c_j)  ;  % number of storage classes
T = size(a_it, 2);  % number of time periods

% start a model
h = rome_model('Linearrule Warehouse Model');

% Define primitive uncertainties
newvar z(M, T) uncertain;
rome_box(z, z_LB, z_UB);

% Define linearrule recourse variables
% num pallets of product i assigned to storage class j in period t (after observing z_(t-1))
v_ijt = rome_empty_var(M, N, T); 
for tt = 1:T
    v_ijt(:, :, tt) = newvar(M, N, z(:, 1:tt-1), 'linearrule');
end

% num pallets of product i retrieved from storage class j in period t (after observing z_t)
w_ijt = rome_empty_var(M, N, T); 
for tt = 1:T
    w_ijt(:, :, tt) = newvar(M, N, z(:, 1:tt), 'linearrule');
end
% num pallets of product i in storage class j at start of period t (after observing t-1)
x_ijt = rome_empty_var(M, N, T+1); 
for tt = 1:T+1
    x_ijt(:, :, tt) = newvar(M, N, z(:, 1:tt-1), 'linearrule');
end

d_it = mean_d_it + z; % uncertain demand vector

% Define Constraints
% Goods Supply Constraint Set
rome_constraint(squeeze(sum(v_ijt, 2)) == a_it);

% Goods Demand Constraint Set
rome_constraint(squeeze(sum(w_ijt, 2)) == d_it);

% Goods Carry Forward Constraint Set
for tt = 1:T
    rome_constraint(x_ijt(:, :, tt+1) == x_ijt(:, :, tt) + v_ijt(:, :, tt) - w_ijt(:, :, tt));
end

% No Initial Inventory Constraint Set
rome_constraint(x_ijt(:, :, 1) == 0);

% No Exceed Capacity Constraint Set
C_jt = repmat(c_j, 1, T);   % expanded matrix of storage capacities
rome_constraint(squeeze(sum(x_ijt(:, :, 1:T) + v_ijt, 1)) <= C_jt);

% Non-negativity constraints
rome_constraint(v_ijt >= 0);
rome_constraint(w_ijt >= 0);
rome_constraint(x_ijt >= 0);

% Cost Objective
S_ijt = repmat(s_j', [M, 1, T]);
R_ijt = repmat(r_j', [M, 1, T]);
g = S_ijt .* strip_rand(v_ijt) + R_ijt .* strip_rand(w_ijt);
g = sum(g(:));
rome_minimize(g);

% Solve
h.solve;

% Output
V = cat(4, h.eval(v_ijt), zeros(M, N, T, T));
W = h.eval(w_ijt);
X = h.eval(x_ijt);
obj_val = h.objective;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
