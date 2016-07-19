% +ROME_EXAMPLES\MODEL_TEST_SEG_INVENTORY Test Script to model and solve
% a simple multi-period inventory model
%
% Tests the effect of using generic segregation.
%
% Modified by:
% 1. Joel (17 Feb 2009)

% display welcome
disp(sprintf('\nTesting Simple Inventory Example (Segregation)...'));

% % define model parameters
% T = 5;  % number of periods
% hold_cost   = 0.02 * ones(T, 1);  % holding cost in each period
% back_cost   = 0.2 * ones(T, 1);  % backorder cost in each period
% back_cost(end) = 10*back_cost(1);
% order_cost  = 0.1 * ones(T, 1); % order cost in each period
% x_MAX       = 260;
% mean_demand = 200;
% z_UB = 20;
% z_LB = -z_UB;
% seg_mean_vec = 15:-1:5;          % segregated mean of primitive uncertainties
% alpha = 0.75;           % parameter for demand process
% beta_vec = 0.8:-0.05:0.2;

% Test parameters provided in DataFile
load('+RomeTest/SegInvData.mat');
inventory_cost = [hold_cost; back_cost]';

% make the LDR patterns for X and Y
pX = logical([tril(ones(T)), zeros(T, 1)]);
pY = logical([ones(T, 1), tril(ones(T))]);

% TEST 0: Changing the segregated support
disp(sprintf('\nTest  0: Sensitivity to segregated mean:'));
x_val_arr = zeros(T, 2*T+1, 2);
res_vec = zeros(1, size(x_val_arr, 3));
for ii = 0:(size(x_val_arr,3) - 1)
    % begin ROME environment
    rome_begin;

    % create a model
    h = rome_model('Segregated Inventory');

    % define variables
    newvar zeta(2*T) uncertain;

    % relate z to zeta
    F = [eye(T), -eye(T)];
    z = F*zeta;

    % specify properties of uncertainties
    if(ii == 0)
        % for first case, specify properties of z
        rome_box(z, z_LB, z_UB);    % support
        z.set_mean(0);              % mean
    else
        % specify properties of zeta
        rome_box(zeta, 0, z_UB);                                % support
        rome_constraint(zeta(1:T) + zeta((T+1):end) <= z_UB);   % support
        zeta.set_mean(10); % mean        
    end

    % define LDRs
    x = rome_linearrule(T, zeta, 'Cone', rome_constants.NNOC, 'Pattern', [pX, pX(:, 2:end)]); % order qty
    y = rome_linearrule(T, zeta, 'Pattern', [pY, pY(:, 2:end)]);  % inventory position

    % define a uncertain demand
    d = (eye(T) + alpha*tril(ones(T), -1)) * z + mean_demand;

    % make the carry-over constraint
    D = eye(T) - diag(ones(T-1,1), -1);
    rome_constraint(D*y - x == -d);

    % make the capacity constraint
    rome_constraint(x <= x_MAX);

%     % make the objective
%     rome_minimize(order_cost' * mean(x) ...
%         + hold_cost' * (rome_supp_bound(y)) ... 
%         + back_cost' * (rome_supp_bound(-y)));
    
    % make the objective
    rome_minimize(order_cost' * mean(x) ...
        + inventory_cost * (rome_supp_bound([y; -y])));

    h.solve; % call solve 
    cur_val = h.objective; % store result    
    x_val_arr(:,:,ii+1) = squeeze(linearpart(h.eval(x)));
    
    curr_err = abs(benchmark_seg_supp_result(ii+1) - cur_val) ...
              + norm(benchmark_seg_supp_result_x(:, :, ii+1) - x_val_arr(:,:,ii+1));
%     curr_err = 0;

    if(ii == 0)
        disp(sprintf('Iter %2d, NoSegSupp, Obj = %0.2f, Err = %4.3f', ...
            ii, cur_val, curr_err));
    else
        disp(sprintf('Iter %2d, SegSupp  , Obj = %0.2f, Err = %4.3f', ...
            ii, cur_val, curr_err));
    end
    res_vec(ii+1) = cur_val;
    rome_end % end environment
end

% TEST 1: Changing the segregated mean
disp(sprintf('\nTest  1: Sensitivity to segregated mean:'));
res_vec = zeros(1, numel(seg_mean_vec));
tic;
for ii = 0:numel(seg_mean_vec)
    if(ii > 0)
        % select parameter
        cur_seg_mean = seg_mean_vec(ii);
    end
    
    % begin ROME environment
    rome_begin;

    % create a model
    h = rome_model('Segregated Inventory');

    % define variables
    newvar zeta(2*T) uncertain;

    % relate z to zeta
    F = [eye(T), -eye(T)];
    z = F*zeta;

    % specify properties of uncertainties
    if(ii == 0)
        % for first case, specify properties of z
        rome_box(z, z_LB, z_UB);        % support
        z.set_mean(0);
    else
        % specify properties of zeta
        rome_box(zeta, 0, z_UB);                    % support
        zeta.set_mean(cur_seg_mean); % mean        
    end

    % define LDRs
    x = rome_linearrule(T, zeta, 'Cone', rome_constants.NNOC, 'Pattern', [pX, pX(:, 2:end)]); % order qty
    y = rome_linearrule(T, zeta, 'Pattern', [pY, pY(:, 2:end)]);  % inventory position

    % define a uncertain demand
    d = (eye(T) + alpha*tril(ones(T), -1)) * z + mean_demand;

    % make the carry-over constraint
    D = eye(T) - diag(ones(T-1,1), -1);
    rome_constraint(D*y - x == -d);

    % make the capacity constraint
    rome_constraint(x <= x_MAX);

    % make the objective
    rome_minimize(order_cost' * mean(x) ...
        + inventory_cost * (rome_supp_bound([y; -y])));

    h.solve; % call solve 
    cur_val = h.objective; % store result    
    rome_end % end environment
    curr_err = abs(benchmark_seg_mean_result(ii+1) - cur_val);
    
    if(ii == 0)
        disp(sprintf('Iter %2d, NoSegMean,      Obj = %0.2f, Err = %4.3f', ...
            ii, cur_val, curr_err));
    else
        disp(sprintf('Iter %2d, SegMean = %4.1f, Obj = %0.2f, Err = %4.3f', ...
            ii, cur_seg_mean, cur_val, curr_err));
    end
    res_vec(ii+1) = cur_val; % store objective value in vector
end
disp(sprintf('Average Time = %0.3f secs', toc / (numel(seg_mean_vec) + 1)));

% TEST 2: Changing the segregated covariance by varying beta
disp(sprintf('\nTest  2: Sensitivity to segregated covariance:'));
res_vec = zeros(1, numel(beta_vec));
tic;
for ii = 0:numel(beta_vec)
    if(ii > 0)
        % select parameter
        beta = beta_vec(ii);
    end
    
    % begin ROME environment
    rome_begin;

    % create a model
    h = rome_model('Segregated Inventory');

    % define variables
    newvar zeta(2*T) uncertain;

    % relate z to zeta
    F = [eye(T), -eye(T)];
    z = F*zeta;

    % specify properties of uncertainties
    if(ii == 0)
        % for first case, specify properties of z
        rome_box(z, z_LB, z_UB);        % support
        z.set_mean(0);   % mean
        z.Covar = (1/3)*(z_UB)^2 * eye(T);   % covar
    else
        % specify properties of zeta
        rome_box(zeta, 0, z_UB);            % support
        zeta.set_mean(10);   % mean
        zeta.Covar = (1/3)*(z_UB)^2 * [beta * eye(T), zeros(T); zeros(T), (1-beta)*eye(T)]; % Covar
    end

    % define LDRs
    x = rome_linearrule(T, zeta, 'Cone', rome_constants.NNOC, 'Pattern', [pX, pX(:, 2:end)]); % order qty
    y = rome_linearrule(T, zeta, 'Pattern', [pY, pY(:, 2:end)]);  % inventory position

    % define a uncertain demand
    d = (eye(T) + alpha*tril(ones(T), -1)) * z + mean_demand;

    % make the carry-over constraint
    D = eye(T) - diag(ones(T-1,1), -1);
    rome_constraint(D*y - x == -d);

    % make the capacity constraint
    rome_constraint(x <= x_MAX);

    % make the objective
    rome_minimize(order_cost' * mean(x) ...
        + inventory_cost * (rome_covar_bound([y; -y])));

    h.solve; % call solve 
    cur_val = h.objective; % store result    
    rome_end % end environment
    curr_err = abs(benchmark_seg_covar_result(ii+1) - cur_val);
    
    if(ii == 0)
        disp(sprintf('Iter %2d, NoSegCovar,     Obj = %0.2f, Err = %4.3f', ...
            ii, cur_val, curr_err));
    else
        disp(sprintf('Iter %2d, SegBeta = %4.2f, Obj = %0.2f, Err = %4.3f', ...
            ii, beta, cur_val, curr_err));
    end
    res_vec(ii+1) = cur_val; % store objective value in vector
end
disp(sprintf('Average Time = %0.3f secs', toc / (numel(beta_vec) + 1)));



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
