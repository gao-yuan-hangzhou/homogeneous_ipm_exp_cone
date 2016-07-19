% ops_test_mtimes.m
% Script to test matrix multiplication

% 1. Testing matrix multiplies
% -----------------------------
disp ('1. Testing Matrix Multiplication ...');

% when variable is a matrix
x = rome_model_var(25, 20);
x_input = rand(size(x));
x_premult = rand(1000, size(x, 1));
x_postmult = rand(size(x, 2), 30);
x_premult_s = rand;
x_postmult_s = rand;

% perform computation
test_x_premult = x_premult * x;
test_x_postmult = x * x_postmult;
test_x_premult_s = x_premult_s * x;
test_x_postmult_s = x * x_postmult_s;
test_x_premult_f = x_premult * full(x);
test_x_postmult_f = full(x) * x_postmult;
test_x_premult = x_premult * x;
test_x_composite_1 = (x_premult * x) * x_postmult;
test_x_composite_2 = x_premult * (x * x_postmult);

% compute errors
x_val = x.insert(x_input);
dx1 = norm(test_x_premult.insert(x_input) - x_premult * x_val);
dx2 = norm(test_x_postmult.insert(x_input) - x_val * x_postmult);
dx1_s = norm(test_x_premult_s.insert(x_input) - x_premult_s * x_val);
dx2_s = norm(test_x_postmult_s.insert(x_input) - x_val * x_postmult_s);
dx1_f = norm(test_x_premult_f.insert(x_input) - x_premult * x_val);
dx2_f = norm(test_x_postmult_f.insert(x_input) - x_val * x_postmult);
dx1_c = norm(test_x_composite_1.insert(x_input) - (x_premult * x_val) * x_postmult);
dx2_c = norm(test_x_composite_2.insert(x_input) - x_premult * (x_val * x_postmult));


% add to output array
test_array = [test_array, dx1, dx2, dx1_s, dx2_s, dx1_f, dx2_f, dx1_c, dx2_c];
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Pre-Mult ', 'Matrix-Matrix Post-Mult', ...
    'Matrix-Scalar Pre-Mult ', 'Matrix-Scalar Post-Mult', ...
    'FMat-Matrix Pre-Mult   ', 'FMat-Matrix Post-Mult  ', ...
    'Composite Matrix 1     ', 'Composite Matrix 2     '};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(20, 1);
    else
        y = rome_model_var(1, 20);
    end    
    y_input = rand(size(y));
    y_premult = rand(1000, size(y, 1));
    y_postmult = rand(size(y, 2), 30);
    y_premult_s = rand;
    y_postmult_s = rand;

    % perform computation
    test_y_premult    = y_premult * y;
    test_y_postmult   = y * y_postmult;
    test_y_premult_s  = y_premult_s * y;
    test_y_postmult_s = y * y_postmult_s;

    y_val = y.insert(y_input);
    dy1 = norm(test_y_premult.insert(y_input) - y_premult * y_val);
    dy2 = norm(test_y_postmult.insert(y_input) - y_val * y_postmult);
    dy1_s = norm(test_y_premult_s.insert(y_input) - y_premult_s * y_val);
    dy2_s = norm(test_y_postmult_s.insert(y_input) - y_val * y_postmult_s);

    % add to output array
    test_array = [test_array, dy1, dy2, dy1_s, dy2_s];
    descr_array = {descr_array{:}, 'Vector-Matrix Pre-Mult ', 'Vector-Matrix Post-Mult', ...
                   'Vector-Scalar Pre-Mult ', 'Vector-Scalar Post-Mult'};
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));
y_premult = rand(1000, 500);
y_postmult = rand(450, 30);
y_premult_v = rand(1000, size(y, 1));
y_postmult_v = rand(size(y, 2), 30);
y_premult_s = rand;
y_postmult_s = rand;

% perform computation
test_y_premult    = y_premult * y;
test_y_postmult   = y * y_postmult;
test_y_premult_v  = y_premult_v * y;
test_y_postmult_v = y * y_postmult_v;
test_y_premult_s  = y_premult_s * y;
test_y_postmult_s = y * y_postmult_s;

y_val = y.insert(y_input);
dy1  = norm(test_y_premult.insert(y_input) - y_premult * y_val);
dy2  = norm(test_y_postmult.insert(y_input) - y_val * y_postmult);
dy1_v = norm(test_y_premult_v.insert(y_input) - y_premult_v * y_val);
dy2_v = norm(test_y_postmult_v.insert(y_input) - y_val * y_postmult_v);
dy1_s = norm(test_y_premult_s.insert(y_input) - y_premult_s * y_val);
dy2_s = norm(test_y_postmult_s.insert(y_input) - y_val * y_postmult_s);

% add to output array
test_array = [test_array, dy1, dy2, dy1_v, dy2_v, dy1_s, dy2_s];
descr_array = {descr_array{:}, ...
    'Scalar-Matrix Pre-Mult ', 'Scalar-Matrix Post-Mult', ...
    'Scalar-Vector Pre-Mult ', 'Scalar-Vector Post-Mult', ...
    'Scalar-Scalar Pre-Mult ', 'Scalar-Scalar Post-Mult'};


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
