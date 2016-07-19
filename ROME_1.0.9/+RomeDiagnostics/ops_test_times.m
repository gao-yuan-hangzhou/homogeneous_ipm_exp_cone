% ops_test_times
% Script To Test Array Multiplication

% 2. Testing array multiplication
% --------------------------------
disp ('2. Testing Array Multiplication ...');

% when variable is a matrix
x = rome_model_var(50, 37);
x_input = rand(size(x));
x_premult = rand(size(x));
x_postmult = rand(size(x));
x_premult_s = rand;
x_postmult_s = rand;

% perform computation
test_x_premult = x_premult .* x;
test_x_postmult = x .* x_postmult;
test_x_premult_s = x_premult_s .* x;
test_x_postmult_s = x .* x_postmult_s;
test_x_premult_f = x_premult .* full(x);
test_x_postmult_f = full(x) .* x_postmult;
test_x_composite_1 = (x_premult .* x) .* x_postmult;
test_x_composite_2 = x_premult .* (x .* x_postmult);

% compute errors
x_val = x.insert(x_input);
dx1 = norm(test_x_premult.insert(x_input) - x_premult .* x_val);
dx2 = norm(test_x_postmult.insert(x_input) - x_val .* x_postmult);
dx1_s = norm(test_x_premult_s.insert(x_input) - x_premult_s .* x_val);
dx2_s = norm(test_x_postmult_s.insert(x_input) - x_val .* x_postmult_s);
dx1_f = norm(test_x_premult_f.insert(x_input) - x_premult .* x_val);
dx2_f = norm(test_x_postmult_f.insert(x_input) - x_val .* x_postmult);
dx1_c = norm(test_x_composite_1.insert(x_input) - (x_premult .* x_val) .* x_postmult);
dx2_c = norm(test_x_composite_2.insert(x_input) - x_premult .* (x_val .* x_postmult));

% add to output array
test_array = [test_array, dx1, dx2, dx1_s, dx2_s, dx1_f, dx2_f, dx1_c, dx2_c];
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Array-Pre-Mult ', 'Matrix-Matrix Array-Post-Mult', ...
    'Matrix-Scalar Array-Pre-Mult ', 'Matrix-Scalar Array-Post-Mult', ...
    'FMat-Matrix Array-Pre-Mult   ', 'FMat-Matrix Array-Post-Mult  ', ...
    'Composite Array Mult 1       ', 'Composite Array Mult 2       '};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(60, 1);
    else
        y = rome_model_var(1, 60);
    end    
    y_input = rand(size(y));
    y_premult = rand(size(y));
    y_postmult = rand(size(y));
    y_premult_s = rand;
    y_postmult_s = rand;

    % perform computation
    test_y_premult    = y_premult .* y;
    test_y_postmult   = y .* y_postmult;
    test_y_premult_s  = y_premult_s .* y;
    test_y_postmult_s = y .* y_postmult_s;

    y_val = y.insert(y_input);
    dy1 = norm(test_y_premult.insert(y_input) - y_premult .* y_val);
    dy2 = norm(test_y_postmult.insert(y_input) - y_val .* y_postmult);
    dy1_s = norm(test_y_premult_s.insert(y_input) - y_premult_s .* y_val);
    dy2_s = norm(test_y_postmult_s.insert(y_input) - y_val .* y_postmult_s);

    % add to output array
    test_array = [test_array, dy1, dy2, dy1_s, dy2_s];
    descr_array = {descr_array{:}, ...
        'Vector-Vector Array-Pre-Mult ', 'Vector-Vector Array-Post-Mult', ...
        'Vector-Scalar Array-Pre-Mult ', 'Vector-Scalar Array-Post-Mult'};
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));
y_premult = rand(700, 550);
y_postmult = rand(300, 30);
y_premult_v = rand(800, 1);
y_postmult_v = rand(1, 340);
y_premult_s = rand;
y_postmult_s = rand;

% perform computation
test_y_premult    = y_premult .* y;
test_y_postmult   = y .* y_postmult;
test_y_premult_v  = y_premult_v .* y;
test_y_postmult_v = y .* y_postmult_v;
test_y_premult_s  = y_premult_s .* y;
test_y_postmult_s = y .* y_postmult_s;

y_val = y.insert(y_input);
dy1  = norm(test_y_premult.insert(y_input) - y_premult .* y_val);
dy2  = norm(test_y_postmult.insert(y_input) - y_val .* y_postmult);
dy1_v = norm(test_y_premult_v.insert(y_input) - y_premult_v .* y_val);
dy2_v = norm(test_y_postmult_v.insert(y_input) - y_val .* y_postmult_v);
dy1_s = norm(test_y_premult_s.insert(y_input) - y_premult_s .* y_val);
dy2_s = norm(test_y_postmult_s.insert(y_input) - y_val .* y_postmult_s);

% add to output array
test_array = [test_array, dy1, dy2, dy1_v, dy2_v, dy1_s, dy2_s];
descr_array = {descr_array{:}, ...
    'Scalar-Matrix Array-Pre-Mult ', 'Scalar-Matrix Array-Post-Mult', ...
    'Scalar-Vector Array-Pre-Mult ', 'Scalar-Vector Array-Post-Mult', ...
    'Scalar-Scalar Array-Pre-Mult ', 'Scalar-Scalar Array-Post-Mult'};



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
