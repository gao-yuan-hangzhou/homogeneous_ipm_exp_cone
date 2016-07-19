% model_test_socp.m
% Simple SOCP to sanity check ROME 

% welcome message
disp(sprintf('\nTesting Simple SOCP Model ... '));

rome_begin;

% SOCP 1: 2-dimensional ball
tic;
h_model_socp = rome_model('First SOCP');
x = rome_model_var;
y = rome_model_var;
rome_maximize(x+y);
rome_constraint(norm2([x;y]) <= 1);
h_model_socp.solve;
x_min = zeros(2, 1);
x_min(1) = h_model_socp.eval(x);
x_min(2) = h_model_socp.eval(y);
t = toc;
disp(sprintf('SOCP 1: Error = %g, time taken = %g secs', norm(x_min(1:2) - [1; 1] ./ sqrt(2)), t));

% SOCP 2: N-dimensional ball
tic;
N = 100;
h_model_socp_2 = rome_model('Second SOCP');
z = rome_model_var(N);
c = rand(size(z));
rome_maximize(c.' * z);
rome_constraint(norm2(z) <= 1);
h_model_socp_2.solve();
x_min = h_model_socp_2.eval(z);
t = toc;
disp(sprintf('SOCP 2: Error = %g, time taken = %g secs', norm(x_min(1:N) - c./norm(c)), t));

% SOCP 3: more complex N-dimensional ball
tic;
N = 10;
h_model_socp_3 = rome_model('Third SOCP');
x = rome_model_var(N);
c = rand(size(x));
L = 2*rand(N)-1;
Q = L * L';

rome_maximize(c.' * x);
rome_constraint(norm2(chol(Q) * x) <= 1);
solve(h_model_socp_3);
x_min = h_model_socp_3.eval(x);
t = toc;

disp(sprintf('SOCP 3: Error = %g, time taken = %g secs', norm(x_min - (inv(Q) * c ./ sqrt((c' * inv(Q) * c) ))), t)); 

% SOCP 4: With constants
tic;
x_sol = 4;
h_model_socp_4 = rome_model('Fourth SOCP');
x = rome_model_var;
rome_maximize(x);
rome_constraint(norm2([x; 3]) <= 5);
solve(h_model_socp_4);
x_max = h_model_socp_4.eval(x);
t = toc;
disp(sprintf('SOCP 4: Error = %g, time taken = %g secs', norm(x_max - x_sol), t)); 

% SOCP 5: With constants
tic;
x_sol = 1;
h_model_socp_5 = rome_model('SOCP with square constraints');
x = rome_model_var(2);
rome_maximize(sum(x));
Q = 4*eye(2);
% Q(1, 2) = 1;
% Q(2, 1) = 1;

rome_constraint(sq(x, Q) <= 2);
h_model_socp_5.solve;
x_max = h_model_socp_5.eval(x);
t = toc;
disp(sprintf('SOCP 5: Error = %g, time taken = %g secs', norm(x_max - x_sol), t)); 



% h_model_socp_2.solve();
% x_min = h_model_socp_2.eval_var(z);
% t = toc;
% disp(sprintf('SOCP 2: Error = %g, time taken = %g secs', norm(x_min(1:N) - c./norm(c)), t));

% % SOCP 3: Test Case
% tic;
% h3 = rome_model('Third SOCP');
% x = rome_model_var(10);
% rome_minimize(sum(x));
% rome_constraint(norm2(x(1:2)) <= 0.001);
% h3.solve;
% t = toc;
% disp(sprintf('SOCP 3: Complete, time taken = %g secs', t));

% rome_end;
% clearvars -except ROME_ENV;



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
