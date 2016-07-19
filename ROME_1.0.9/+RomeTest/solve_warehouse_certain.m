function [obj_val, V, W, X] = solve_warehouse_certain(z, mean_d_it, a_it, c_j, s_j, r_j)

% 
% +ROMETEST\SOLVE_WAREHOUSE_CERTAIN Helper routine to solve
% deterministic warehouse problem.
%
%   [obj_val, V, W, X] = SOLVE_WAREHOUSE(z, d, A, C, S, R) returns the
%   minimum cost, storage policy V, retrieval policy W, and holding policy
%   X for the supplied data.
%   
%   Define: M as the number of types of good, N as number of storage
%   classes, and T as the planning horizon, the arguments to this function
%   are:
%
%   z: M by T matrix of instantiated primitive uncertainties.
%   d: M by T matrix of mean demand in each period. The actual demand is
%      computed as d_actual = d + z;
%   A: M by T matrix of fixed supply of goods in each period.
%   C: N by 1 column vector of capacity of each storage class
%   S: N by 1 column vector of storage costs for the storage classes
%   R: N by 1 column vector of retrieval costs for the storage classes
%
%   Defaults: M = 3, T = 3, N = 5. Default values for each of the arguments
%   can be found in the .m file.
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

% begin rome environment
rome_begin;

% default arguments
if(nargin < 6) 
    r_j  = [  1, 1.5,  2,  3, 100]';      % column-vector of retrieve costs
end
if(nargin < 5) 
    s_j  = [  1, 1.5,  2,  3, 100]';      % column-vector of store costs
end
if(nargin < 4) 
    c_j  = [ 15,  30, 45, 60, 300000]';   % column-vector of storage capacities
end
if(nargin < 3) 
    a_it = [35 35 35; ...   % fixed arrival of product i at time t
            35 35 35; ...
            10 10 25];  
end
if(nargin < 2) 
    mean_d_it = [29, 19, 45; ... % mean demand of product i at time t
                 10, 47, 37; ...
                  6,  7, 20];
end
if(nargin < 1) 
    z = zeros(size(a_it));            % z represents an instant of a uncertainty
end

% DATA
M = size(a_it, 1);  % number of products
N = length(c_j)  ;  % number of storage classes
T = size(a_it, 2);  % number of time periods

% start a model
h = rome_model('Deterministic Warehouse Model');

% define demand
d_it = mean_d_it + z;

% Define decision variables
v_ijt = rome_model_var(M, N, T, 'Cone', rome_constants.NNOC);  % num pallets of product i assigned to storage class j in period t
w_ijt = rome_model_var(M, N, T, 'Cone', rome_constants.NNOC);  % num pallets of product i retrieved from storage class j in period t
x_ijt = rome_model_var(M, N, T+1 , 'Cone', rome_constants.NNOC); % num pallets of product i in storage class j at start of period t

% Define Constraints
% ------------------
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
C_jt = repmat(c_j, 1, T);
rome_constraint(squeeze(sum(x_ijt(:, :, 1:T) + v_ijt, 1)) <= C_jt);

% Define Cost objective
% ------------------
S_ijt = repmat(s_j', [M, 1, T]);        % expand cost vectors into 3-D arrays
R_ijt = repmat(r_j', [M, 1, T]);
g = S_ijt .* v_ijt + R_ijt .* w_ijt;    % array multiply and sum
g = sum(g(:));
rome_minimize(g);                       % input into objective

% Solve
% ------
h.solve;

% put values into output
V = h.eval(v_ijt);
W = h.eval(w_ijt);
X = h.eval(x_ijt);
obj_val = h.objective;


% Terminate rome environment
rome_end;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
