% ROMETEST\model_test_norm_min.m

% script for testing some old rome functions
load +RomeTest/NormMinData;

% Welcome Message
disp(sprintf('\nTesting Norm Minimization ... '));

rome_begin;

tic;
h = rome_model('MinNorm');
x_opt = rome_model_var(length(x));
rome_minimize(norm2(x_opt));
rome_constraint(A*x_opt == b);
rome_constraint((norm2(x_opt(1:2:100))) + norm2(x_opt(2:2:100)) <= 0.0001);
h.solve;
x_rome = h.eval(x_opt);
t = toc;

% display results
disp(sprintf('Error vs old = %g, Time Taken = %g', norm(x_rome - x_old_prof), t));

rome_end;

% clear vars
clearvars -except ROME_ENV;

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
