% +ROME_TEST\MODEL_TEST_WAREHOUSE Test Script to model and solve
% the warehouse problem
%
% Calls the following subroutines to solve the actual problem:
% 1. SOLVE_WAREHOUSE_CERTAIN (solves the deterministic problem under
%    perfect information).
% 2. SOLVE_WAREHOUSE_LINEARRULE (solves the robust problem using a
%    linear-decision rule)
%
% References:
% 1. 
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

import RomeTest.*;

% display welcome
disp(sprintf('\nTesting Warehouse Example...'));

% define model parameters
R  = [  1, 1.5,  2,  3, 100]';      % column-vector of retrieve costs
S  = [  1, 1.5,  2,  3, 100]';      % column-vector of store costs
C  = [ 15,  30, 45, 60, 300000]';   % column-vector of storage capacities
A  = [35 35 35; ...     % fixed arrival of product i at time t
      35 35 35; ...
      10 10 25];  
D  = [29, 19, 45; ...    % mean demand of product i at time t
      10, 47, 37; ...
       6,  7, 20];

M = size(A, 1);  % number of products
N = length(C) ;  % number of storage classes
T = size(A, 2);  % number of time periods

% primitive uncertainty representation
zLB = -1;
zUB = 1;

% 1. Solve using linearrule
tic; [E_cost, V_ldr, W_ldr, X_ldr] = solve_warehouse_linearrule(zLB, zUB, D, A, C, S, R); t = toc;
disp('1. Solving using LDR... ');
disp(sprintf('E(Obj) = %0.2f, Time = %0.4f secs\n', E_cost, t));

% 2. Solve deterministic (perfect info model)
SS = repmat(S', [M, 1, T]); % expanded storage costs
RR = repmat(R', [M, 1, T]); % expanded retrieve costs

% define interation parameters
% A) random testing
Niter1 = 20;
all_z1 = unifrnd(zLB, zUB, [size(A), Niter1]);

% B) Extrema testing
% Niter2 = 2^(M*T);
% all_z2 = zeros([size(A), Niter2]);
% for ii = 1:Niter2
%     z = sscanf(dec2bin(ii-1), '%1d');
%     z = [zeros(M*T - numel(z), 1); z];
%     all_z2(:,:,ii) = reshape(z, M, T);
% end
% all_z2 = (zUB - zLB) * all_z2 + zLB;

% Combine all z's together
% all_z = cat(3, all_z1, all_z2);
% Niter = Niter1 + Niter2

% HACK: Just test random only
Niter = Niter1;
all_z = all_z1;

disp('2. Solving under perfect information ... ');
for ii = 1:Niter
    % instantiate primitive uncertainty
    z = all_z(:, :, ii);
    
    % call solving subroutine
    tic; [obj_val, V, W, X] = solve_warehouse_certain(z, D, A, C, S, R); t = toc;

    % compare with linearrule
    v_predict = reshape(reshape(V_ldr, [M*N*T, 1 + M*T])     * [1; z(:)], M, N, T);
    w_predict = reshape(reshape(W_ldr, [M*N*T, 1 + M*T])     * [1; z(:)], M, N, T);
    x_predict = reshape(reshape(X_ldr, [M*N*(T+1), 1 + M*T]) * [1; z(:)], M, N, T + 1);
    d_predict = D + z;
    
    % feasibility check
    feas_vec = check_warehouse_constraints(v_predict, w_predict, x_predict, d_predict, A, C);
    
    % evaluate linearrule cost
    ldr_cost = SS .* v_predict + RR .* w_predict;
    ldr_cost = sum(ldr_cost(:));
    
    if(all(feas_vec))
        % make report    
        disp(sprintf('Iter %3d, PI Obj = %0.2f, LDR Obj = %0.2f, Eff = %0.4f%% Time = %0.4f secs', ...
            ii, obj_val, ldr_cost, obj_val ./ ldr_cost * 100 , t));
    else
        % make error report
        error('Iter %3d FAILED, PI Obj = %0.2f, LDR Obj = %0.2f, Eff = %0.4f%% Time = %0.4f secs', ...
            ii, obj_val, ldr_cost, obj_val ./ ldr_cost * 100 , t);
    end
end

% clear the workspace
clearvars -except ROME_ENV


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
