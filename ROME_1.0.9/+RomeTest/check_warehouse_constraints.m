function feas_vec = check_warehouse_constraints(v, w, x, d, A, C)

% 
% +PROFEXAMPLES\CHECK_WAREHOUSE_CONSTRAINTS Helper routine to check
% feasibility of a particular solution
%
%   feas_vec = SOLVE_WAREHOUSE(v, w, x d, A, C) returns a logical 8x1
%   vector which represents feasibility of different sets of inequalities
%   in the warehouse problem. See the .m file for which inequality each
%   index corresponds to.
%   
%   Define: M as the number of types of good, N as number of storage
%   classes, and T as the planning horizon, the arguments to this function
%   are:
%
%   v: M by N by T matrix representing storage policy of each good to each
%      storage class in each time period
%   w: M by N by T matrix representing retrieval policy of each good from
%      each storage class in each time period
%   d: M by T matrix of actual demand for each good in each time period
%   A: M by T matrix of fixed supply of goods in each period.
%   C: N by 1 column vector of capacity of each storage class
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

M = size(A, 1);  % number of products
N = length(C) ;  % number of storage classes
T = size(A, 2);  % number of time periods
feas_vec = true(8, 1);

CC = repmat(C, 1, T); % expanded capacities

% holdover constraints
for tt = 1:T
    L = x(:,:, tt+1) - (x(:,:,tt) + v(:,:,tt) - w(:,:,tt)) <= 1E-8;
    feas_vec(1) = feas_vec(1) && all(L(:));
end

% no initial inventory constraint
L = (abs(x(:, :, 1)) < 1E-8);
feas_vec(2) = all(L(:));

% Supply constraint
L = abs(squeeze(sum(v, 2)) - A) <= 1E-8;
feas_vec(3) = all(L(:));

% Demand constraint
L = abs(squeeze(sum(w, 2)) - d) <= 1E-8;
feas_vec(4) = all(L(:));

% Capacity constraint
L = (squeeze(sum(x(:,:, 1:T) + v, 1)) - CC) <= 1E-8;
feas_vec(5) = all(L(:));

% Nonnegativity Constraints
L = (v >= 0);
feas_vec(6) = all(L(:));

L = (w >= -1E-8);
feas_vec(7) = all(L(:));

L = (x >= 0);
feas_vec(8) = all(L(:));


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
