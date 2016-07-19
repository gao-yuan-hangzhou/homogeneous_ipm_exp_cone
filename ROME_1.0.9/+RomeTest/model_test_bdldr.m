% model_test_bdldr
% parameters
N = 1;

% Welcome Message
disp(sprintf('\nBegin BDLDR test for N = %d ... ', N));

disp(sprintf('\nUsing LDR'));
h = rome_begin('LDR Base Model');
tic;
% uncertainties
z = rome_rand_model_var(N);
z.set_mean(0);
z.Covar = eye(N);

% model variables
u = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
v = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
s = rome_linearrule(1, z);

% constraints
rome_constraint(u - v == s - sum(z));
rome_box(s, 0, 1);

% objective
rome_minimize(mean(u+v));

% solve
h.solve; 
disp(sprintf('LDR Obj = %0.4f, time = %0.2f secs', h.objective, toc));

% DLDR 
% -----
disp(sprintf('\nUsing DLDR'));
h = rome_begin('DLDR Base Model');
tic;

% uncertainties
z = rome_rand_model_var(N);
z.set_mean(0);
z.Covar = eye(N);

% model variables
u = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
v = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
s = rome_linearrule(1, z);

% constraints
rome_constraint(u - v == s - sum(z));
rome_box(s, 0, 1);

% objective
rome_minimize( mean(u+v) );

% to be placed after objective
[X_vals, f_coeffs] = h.apply_dldr(u + v, [u;v], @rome_covar_bound);
% [X_vals, f_coeffs] = h.apply_na_bdldr(u + v, [u;v], @rome_covar_bound);

% solve
h.solve;
disp(sprintf('DLDR Obj = %0.4f, Err = %0.4f, time = %0.2f secs', ...
     h.objective, abs(h.objective - sqrt(N)), toc));
 
% [deflect_coeffs, deflect_values] =  h.eval_deflect(u)
% [deflect_coeffs, deflect_values] =  h.eval_deflect(v)
% [deflect_coeffs, deflect_values] =  h.eval_deflect(s)
 
% s = squeeze(h.eval_var(s))'; display(s);
% u = squeeze(h.eval_var(u))'; display(u);
% v = squeeze(h.eval_var(v))'; display(v);

% Direct DLDR 
% -----------
disp(sprintf('\nUsing Direct DLDR'));
h = rome_begin('Direct DLDR Base Model');
tic;

% deterministic
u = rome_model_var(N+1);
v = rome_model_var(N+1);
s0 = rome_model_var(1);

% constraints
rome_constraint(u(1) - v(1) == s0);
rome_constraint(u(2:end) - v(2:end) == -1);
rome_box(s0, 0, 1);

% expectation constraint
rome_minimize(norm2(u) + norm2(v));

% solve
h.solve;
disp(sprintf('DLDR Obj = %0.4f, Err = %0.4f, time = %0.2f secs', ...
     h.objective, abs(h.objective - sqrt(N)), toc));
% s0 = squeeze(h.eval_var(s0))'; display(s0);
% u = squeeze(h.eval_var(u))'; display(u);
% v = squeeze(h.eval_var(v))'; display(v);

% BDLDR 
% -----
disp(sprintf('\nUsing BDLDR'));
h = rome_begin('BDLDR Base Model');
tic;

% uncertainties
z = rome_rand_model_var(N);
z.set_mean(0);
z.Covar = eye(N);

% model variables
u = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
v = rome_linearrule(1, z, 'Cone', rome_constants.NNOC);
s = rome_linearrule(1, z);
rome_box(s, 0, 1);

% constraints
rome_constraint(u - v == s - sum(z));

% expectation constraint
rome_minimize(mean(u+v));

% solve
h.solve_deflected;
disp(sprintf('BDLDR Obj = %0.4f, Err = %0.4f, time = %0.2f secs', ...
     h.objective, abs(h.objective - 0.5*(sqrt(N) + sqrt(N+1) - 1)), toc));

% [deflect_coeffs, deflect_vals] = h.eval_deflect(u)
% [deflect_coeffs, deflect_vals] = h.eval_deflect(v)
% [deflect_coeffs, deflect_vals] = h.eval_deflect(s)
 
% s = squeeze(h.eval_var(s))'; display(s);
% u = squeeze(h.eval_var(u))'; display(u);
% v = squeeze(h.eval_var(v))'; display(v);
 
% Direct BDLDR 
% -----------
disp(sprintf('\nUsing Direct BDLDR'));
h = rome_begin('direct BDLDR Base Model');
tic;

% deterministic
u = rome_model_var(N+1);
v = rome_model_var(N+1);

% constraints
s1 = u - v + [0; ones(N, 1)];
s2 = u - v + [-1; ones(N, 1)];

% objective
rome_minimize(norm2(u) + norm2(v) - 0.5 + 0.5*norm2(s1) + 0.5*norm2(s2));

% solve
h.solve;
disp(sprintf('BDLDR Obj = %0.4f, Err = %0.4f, time = %0.2f secs', ...
     h.objective, abs(h.objective - 0.5*(sqrt(N) + sqrt(N+1) - 1)), toc));
% u = squeeze(h.eval_var(u))'; display(u);
% v = squeeze(h.eval_var(v))'; display(v);
% s1 = squeeze(h.eval_var(s1))'; display(s1);
% s2 = squeeze(h.eval_var(s2))'; display(s2);

% Dual solution BDLDR 
% --------------------
disp(sprintf('\nUsing Dual Solution (N = 1)'));
h = rome_begin('Dual BDLDR');
tic;

% deterministic
s = rome_model_var(3);

% objective
rome_minimize(s(1) + s(3));

% SOC constraints in hypoquad
hypoquad(s(1)    , s(3)       , 0.5*(s(2) + 1));
hypoquad(s(1)    , s(3)       , 0.5*(s(2)));
hypoquad(s(1) + 1, s(3)       , 0.5*(s(2) - 1));

% solve
h.solve;
disp(sprintf('Dual BDLDR Obj = %0.4f, time = %0.2f secs', h.objective, toc));

% Correct solution
correct_objective = h.objective;


% Using a tighter bound
% ------------------------
disp(sprintf('\nUsing BDLDR with Tighter Bound'));
h = rome_begin('BDLDR Tighter Bound');
tic;

% uncertainties
z = rome_rand_model_var(N);
z.set_mean(0);
z.Covar = eye(N);

% model variables
u = rome_linearrule(1, z);
v = rome_linearrule(1, z);
s = rome_linearrule(1, z);

% constraints
rome_constraint(u - v == s - sum(z));

% expectation constraint
rome_minimize(sum(rome_pwl_create_bound([u; v], [0,0], [-1, 1])) + ...
                rome_pwl_create_bound(s, [-1, 0, 0], [1, 0, -1]));
              

% rome_minimize(rome_covar_bound(u) + rome_covar_bound(-u) +... 
%               rome_covar_bound(v) + rome_covar_bound(-v) +...               
%               rome_pwl_covar_bound(s, [0, -1; 0, 0]) + ...
%               rome_covar_bound(s-1));
          
% solve
h.solve;
disp(sprintf('BDLDR Obj = %0.4f, Err = %0.4f, time = %0.2f secs', ...
     h.objective, abs(h.objective - correct_objective), toc));


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
