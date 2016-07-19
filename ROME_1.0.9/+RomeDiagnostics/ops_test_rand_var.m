% ops_test_rand_var.m
% Script to test matrix multiplication for uncertain variables

% 8. Testing matrix multiplies
% -----------------------------
disp ('8. Testing Uncertain Variable Matrix Multiplication ...');

% parameters 
N = 20; M = 30;

% model variables
x = rome_model_var(N, M, 'Cone', rome_constants.NNOC);
y = rome_model_var(N, M, 'Cone', rome_constants.NNOC);
p = rome_rand_model_var(N, M);
q = rome_rand_model_var(N, M);

% define input values
x_input = rand(N, M);
y_input = rand(N, M);
p_input = rand(N, M);
q_input = rand(N, M);

% INSERTION
% insert values
dp1 = norm(p.insert(p_input) - p_input, 'fro');
dq1 = norm(q.insert(q_input) - q_input, 'fro');
dx1 = norm(x.insert(x_input) - x_input, 'fro');
dy1 = norm(y.insert(y_input) - y_input, 'fro');

test_array = [test_array, dp1, dq1, dx1, dy1];
descr_array = {descr_array{:}, ...
    'Pure Rand Insert1   ', 'Pure Rand Insert2   ', ...
    'Pure Certain Insert1', 'Pure Certain Insert2', ...
    };

% INNER PRODUCT
% Multiply 
r1 = p.' * x;
r2 = x.' * p;

s1 = q.' * y;
s2 = y.' * q;

% Multiply (exchange terms)
t1 = q.' * x;
t2 = x.' * q;

u1 = p.' * y;
u2 = y.' * p;

% compute true values
r1_val = p_input.' * x_input;
r2_val = x_input.' * p_input;

s1_val = q_input.' * y_input;
s2_val = y_input.' * q_input;

t1_val = q_input.' * x_input;
t2_val = x_input.' * q_input;

u1_val = p_input.' * y_input;
u2_val = y_input.' * p_input;

% compute errors
dr1 = norm(r1.insert(x_input, p_input) - r1_val, 'fro');
dr2 = norm(r2.insert(x_input, p_input) - r2_val, 'fro');

ds1 = norm(s1.insert(y_input, q_input) - s1_val, 'fro');
ds2 = norm(s2.insert(y_input, q_input) - s2_val, 'fro');

dt1 = norm(t1.insert(x_input, q_input) - t1_val, 'fro');
dt2 = norm(t2.insert(x_input, q_input) - t2_val, 'fro');

du1 = norm(u1.insert(y_input, p_input) - u1_val, 'fro');
du2 = norm(u2.insert(y_input, p_input) - u2_val, 'fro');

% Append to Output Array
test_array = [test_array, dr1, dr2, ds1, ds2, dt1, dt2, du1, du2];
descr_array = {descr_array{:}, ...
    'Rand*Certain Inner Insert1', 'Certain*Rand Inner Insert1', ...
    'Rand*Certain Inner Insert2', 'Certain*Rand Inner Insert2', ...
    'Rand*Certain Inner Insert3', 'Certain*Rand Inner Insert3', ...
    'Rand*Certain Inner Insert4', 'Certain*Rand Inner Insert4', ...
    };

% OUTER-PRODUCT
% Multiply
r1 = p * x.';
r2 = x * p.';

s1 = q * y.';
s2 = y * q.';

% Multiply (inner product, exchange terms)
t1 = q * x.';
t2 = x * q.';

u1 = p * y.';
u2 = y * p.';

% compute true values
r1_val = p_input * x_input.';
r2_val = x_input * p_input.';

s1_val = q_input * y_input.';
s2_val = y_input * q_input.';

t1_val = q_input * x_input.';
t2_val = x_input * q_input.';

u1_val = p_input * y_input.';
u2_val = y_input * p_input.';


% compute errors
dr1 = norm(r1.insert(x_input, p_input) - r1_val, 'fro');
dr2 = norm(r2.insert(x_input, p_input) - r2_val, 'fro');

ds1 = norm(s1.insert(y_input, q_input) - s1_val, 'fro');
ds2 = norm(s2.insert(y_input, q_input) - s2_val, 'fro');

dt1 = norm(t1.insert(x_input, q_input) - t1_val, 'fro');
dt2 = norm(t2.insert(x_input, q_input) - t2_val, 'fro');

du1 = norm(u1.insert(y_input, p_input) - u1_val, 'fro');
du2 = norm(u2.insert(y_input, p_input) - u2_val, 'fro');

% Append to Output Array
test_array = [test_array, dr1, dr2, ds1, ds2, dt1, dt2, du1, du2];
descr_array = {descr_array{:}, ...
    'Rand*Certain Outer Insert1', 'Certain*Rand Outer Insert1', ...
    'Rand*Certain Outer Insert2', 'Certain*Rand Outer Insert2', ...
    'Rand*Certain Outer Insert3', 'Certain*Rand Outer Insert3', ...
    'Rand*Certain Outer Insert4', 'Certain*Rand Outer Insert4', ...
    };


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
