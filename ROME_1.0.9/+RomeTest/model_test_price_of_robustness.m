% +ROMETEST\MODEL_TEST_PRICE_OF_ROBUSTNESS Test Script to model and solve
% a simple portfolio problem. 
%
% This uses three different methods of solving the problem
% 1. L1 - norm method
% 2. Describing Beta (defined in the Paper) by a linear decision rule
% 3. Direct translation of the dual (fully deterministic formulation)
%
% References:
% 1. D. Bertsimas, M.Sim. "The Price of Robustness", Operations Research,
%    2003, 54(1), 35-53. 
%
% Modified by:
% 1. Joel (22 Oct 2008)
%

% Begin Rome Model
disp(sprintf('\nTesting Price of Robustness ... '));

% Parameters
N = 150;
ind = (1:N)';
avg_p = 1.15 + ind * (0.05 / 150);
sigma = (0.05/450) * sqrt(2 * ind * N * (N + 1));

% Range of values of gamma
Gamma = 0:5:150;
robust_return_L1   = zeros(numel(Gamma), 1);
robust_return_beta = zeros(numel(Gamma), 1);
robust_return_dual = zeros(numel(Gamma), 1);

exp_return_L1   = zeros(numel(Gamma), 1);
exp_return_beta = zeros(numel(Gamma), 1);
exp_return_dual = zeros(numel(Gamma), 1);

% iterate over each Gamma
disp(sprintf('\nUsing L1 Norm ...'));
for ii = 1:numel(Gamma)
    % Start Model
    rome_begin;
    h = rome_model('Price of Robustness (using L1)');
    
    % display
    disp(sprintf('Iter %3d, Gamma = %0.1f ', ii, Gamma(ii)));
    
    % Portfolio weights (don't allow shorting)
    x = rome_model_var(N, 'Cone', rome_constants.NNOC);
    rome_constraint(sum(x) == 1);
    
    % uncertain portfolio returns
    p = rome_rand_model_var(N);
    rome_constraint(p <= avg_p + sigma);
    rome_constraint(p >= avg_p - sigma);
    rome_constraint(norm1((p - avg_p) ./ sigma) <= Gamma(ii));
  
    % robust portfolio return
    r = rome_model_var;
    rome_constraint(r <= p' * x );
    rome_maximize(r);

    % solve
    h.solve;

    x_val = h.eval(x);
    exp_return_L1(ii) = avg_p' * x_val;
    robust_return_L1(ii) = h.objective;
    
    % end model
    rome_end;
end

% iterate over each Gamma
disp(sprintf('\nUsing Beta ...'));
for ii = 1:numel(Gamma)
    % Start Model
    rome_begin;
    h = rome_model('Price of Robustness (using Beta)');
    
    % display
    disp(sprintf('Iter %3d, Gamma = %0.1f ', ii, Gamma(ii)));
    
    % Portfolio weights (don't allow shorting)
    x = rome_model_var(N, 'Cone', rome_constants.NNOC);
    rome_constraint(sum(x) == 1);
    
    % uncertain portfolio returns
    z = rome_rand_model_var(N);
    rome_constraint(sum(z) <= Gamma(ii));
    rome_constraint(z >= 0);
    rome_constraint(z <= 1);    
    
    % robust portfolio return
    r = rome_model_var;
    rome_constraint(r <= (avg_p' * x) - sum(sigma .* x .* z));
    rome_maximize(r);

    % solve
    h.solve;

    x_val = h.eval(x);
    exp_return_beta(ii) = avg_p' * x_val;
    robust_return_beta(ii) = h.objective;
    
    % end model
    rome_end;
end

% iterate over each Gamma
disp(sprintf('\nUsing Direct Dual ...'));
for ii = 1:numel(Gamma)
    % Start Model
    rome_begin;
    h = rome_model('Price of Robustness (explicit dual)');
    
    % display
    disp(sprintf('Iter %3d, Gamma = %0.1f ', ii, Gamma(ii)));
    
    % Portfolio weights (don't allow shorting)
    x = rome_model_var(N, 'Cone', rome_constants.NNOC);
    rome_constraint(sum(x) == 1);
    
    % variables from duality
    y = rome_model_var(N, 'Cone', rome_constants.NNOC);
    z = rome_model_var('Cone', rome_constants.NNOC);
    rome_constraint(y + z >= sigma .* x);
    
    % robust return
    r = rome_model_var('Cone', rome_constants.NNOC);
    rome_constraint(r <= avg_p' * x - (z * Gamma(ii) + sum(y)));
    rome_maximize(r);
    
    % solve
    h.solve;

    x_val = h.eval(x);
    exp_return_dual(ii) = avg_p' * x_val;
    robust_return_dual(ii) = h.objective;
    
    % end model
    rome_end;
end

load +RomeTest/PriceRobustData.mat

% Error Check:
robust_diff = bsxfun(@minus, [robust_return_L1, robust_return_beta, robust_return_dual], robust_return_benchmark);
exp_diff = bsxfun(@minus, [exp_return_L1, exp_return_beta, exp_return_dual], exp_return_benchmark);

robust_err = sum(robust_diff.^2, 2);
exp_err = sum(exp_diff.^2, 2);

% report
disp(sprintf('L1   method, robust err = %0.2f, expected err = %0.2f', robust_err(1), exp_err(1)));
disp(sprintf('Beta method, robust err = %0.2f, expected err = %0.2f', robust_err(2), exp_err(2)));
disp(sprintf('Dual method, robust err = %0.2f, expected err = %0.2f', robust_err(3), exp_err(3)));

% remove all unnecessary variables
clearvars -except ROME_ENV


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
