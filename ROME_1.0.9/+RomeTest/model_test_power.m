% model_test_power.m
% Test Script to Check Power and Exp functionality in ROME

% % check here
% h = rome_begin('test');
% x = rome_model_var(10);
% y = powercone(x, [1, 2], 2);
% return;

% 1. Testing Powers
% --------------
% define range of n
n_range = 1:2:30;
x_vec = zeros(numel(n_range), 1);
err_vec = zeros(numel(n_range), 1);

% Display welcome message
disp(sprintf('\nTesting Powers from %d to %d, skip = %d ...', ...
    n_range(1), n_range(end), n_range(2) - n_range(1)));

tic;
for ii = 1:numel(n_range)
    % begin
    rome_begin;

    % make model
    h = rome_model('Test Powers');

    % create variable
    x = rome_model_var('Cone', rome_constants.NNOC);

    % make power consrtaint
    rome_constraint(x.^n_range(ii) <= 2);
    rome_maximize(x);

    h.solve;
    x_vec(ii) = h.eval(x);
    err_vec(ii) = norm(2^(1/n_range(ii)) - x_vec(ii));

    % end
    rome_end;

end
t= toc;

% display results
disp(sprintf('Power Test Complete, Net Error = %0.3f. Average time taken = %0.2f secs', sum(err_vec) ./ numel(n_range), t));

% 1B. Testing Powers (B)
% ------------------------
b_vec   = [5, -5];
N = length(b_vec);
err_vec = zeros(N, 1);
for ii = 1:N
    h = rome_begin;
    newvar x;
    rome_minimize((x - b_vec(ii))^2); 
    h.solve;
    err_vec = h.eval(x) - b_vec(ii);
end
disp(sprintf('Power Test B Complete, Net Error = %0.3f. Average time taken = %0.2f secs', sum(err_vec) ./ N, t));

% Testing Exponential 
% ---------------------
disp(sprintf('\nTesting Exponential ...'));
rome_begin;

h = rome_model('Test Exponential');
tic;
x = rome_model_var('Cone', rome_constants.NNOC);
rome_constraint(exp(x) <= 2);
rome_maximize(x);
h.solve;
t = toc;

% check
cur_err = norm(log(2) - h.objective);
disp(sprintf('Exp Test Complete. Net Error = %0.3f, Average Time taken = %0.2f secs', cur_err, t));

% complete
rome_end;

% Testing Power Cone
% ------------------------
disp(sprintf('\nTesting Power Cone ...'));
rome_begin;

% Parameters
n_range = 2:3:45;
N = numel(n_range);
C = 15;
M = 10;
cone_bound = unidrnd(C, N, 1) + C / 2;
mu_bound   = unidrnd(M, N, 1);
err_vec    = zeros(N, 1);

% Iterate
tic;
for ii = 1:N
    h = rome_model('Test Power Cone');

    x  = rome_model_var('Cone', rome_constants.NNOC);
    mu = rome_model_var('Cone', rome_constants.NNOC);
    rome_constraint(mu <= mu_bound(ii));
    rome_constraint(powercone(x, mu, n_range(ii)) <= cone_bound(ii));
    rome_maximize(x);
    
    % Solve
    h.solve;
    x_val = h.eval(x);
    mu_val = h.eval(mu);

    % check error
    err_vec(ii) = pos(mu_val * (x_val / mu_val).^n_range(ii) - cone_bound(ii)) + pos(mu_val - mu_bound(ii));
end
t = toc;

% check
disp(sprintf('Power Cone Test Complete. Net Error = %0.6f, Average Time taken = %0.2f secs', norm(err_vec), t / N));
rome_end;

% Testing Exponential Cone
% ------------------------
disp(sprintf('\nTesting Exponential Cone Part 1...'));
rome_begin;

% Parameters
N = 20;
C = 15;
M = 10;
cone_bound = unidrnd(C, N, 1) + C / 2;
mu_bound   = unidrnd(M, N, 1);
err_vec    = zeros(N, 1);

% Iterate
tic;
for ii = 1:N
    h = rome_model('Test Exponential Cone');

    x  = rome_model_var;
    mu = rome_model_var('Cone', rome_constants.NNOC);
    rome_constraint(mu <= mu_bound(ii));
    rome_constraint(expcone(x, mu) <= cone_bound(ii));
    rome_maximize(x);
    
    % Solve
    h.solve;
    x_val = h.eval(x);
    mu_val = h.eval(mu);

    % check error
    err_vec(ii) = pos(mu_val * exp(x_val / mu_val) - cone_bound(ii)) + pos(mu_val - mu_bound(ii));
end
t = toc;

% check
disp(sprintf('Exp Cone Test Part 1 Complete. Net Error = %0.6f, Average Time taken = %0.2f secs', norm(err_vec), t / N));
rome_end;

% Testing Exponential Cone PART 2
% -------------------------------
disp(sprintf('\nTesting Exponential Cone Part 2...'));
rome_begin;

% Parameters
N = 20;
a_arr = linspace(-20,20, N);
true_arr = zeros(N, 1);
opt_arr  = zeros(N, 1);
b_val = 1;

for ii = 1:N
    a = a_arr(ii);
    h = rome_model('Test Exponential Cone');

    % make variables
    x  = rome_model_var;
    mu = rome_model_var('Cone', rome_constants.NNOC);
    d  = rome_model_var;
    b = hypoquad(mu, d);
    c  = rome_model_var;
    a_var = rome_model_var;

    % make constraints
    rome_constraint(a_var == a);
    rome_constraint(b == b_val);
    rome_constraint(a_var + d <= x);
    rome_constraint(expcone(x, mu) <= c);

%test quad
%     mu = rome_model_var('Cone', rome_constants.NNOC);
%     a_var = rome_model_var;
%     b = rome_model_var;
%     c = rome_model_var;
%     rome_constraint(b == b_val);
%     rome_constraint(a_var == a);
%     rome_constraint(quadexpcone(a, b, mu) <= c);

    
    % minimize epigraph bound
    rome_minimize(c);

    % Solve
    h.solve;

    % get values
    x_val = h.eval(x);
    mu_val = h.eval(mu);
%     d_val  = h.eval(d);
%     b_val  = h.eval(b);
%     c_val  = h.eval(c); 

    % check error
    true_mu = (a + sqrt(a.^2 + 8*b_val.^2)) / 2;
    true_val = true_mu * exp(a./true_mu + 1./true_mu.^2);
    true_arr(ii) = true_val;

    mu_err = true_mu - mu_val;
    err = true_val - h.objective;
    opt_arr(ii) = h.objective;
end

% report:
figure;
semilogy(a_arr, opt_arr, 'r-', a_arr, true_arr, 'b-');
% legend('Approx. Value', 'Exact Value');
disp(sprintf('Exp Cone Test Part 2 Complete'));


% clear vals 
rome_end;
clearvars -except ROME_ENV;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
