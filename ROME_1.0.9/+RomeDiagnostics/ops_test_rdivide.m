% ops_test_rdivide
% Script To Test Array Division

% 3. Testing array division
% --------------------------------
disp ('3. Testing Array Division ...');

% when variable is a matrix
x = rome_model_var(40, 32);
x_input = rand(size(x));
x_postdiv = rand(size(x));
x_postdiv_s = rand;

% perform computation
test_x_postdiv = x ./ x_postdiv;
test_x_postdiv_s = x ./ x_postdiv_s;
test_x_postdiv_f = full(x) ./ x_postdiv;
test_x_postdiv_fs = full(x) ./ x_postdiv_s;

% compute errors
x_val = x.insert(x_input);
dx2 = norm(test_x_postdiv.insert(x_input) - x_val ./ x_postdiv);
dx2_s = norm(test_x_postdiv_s.insert(x_input) - x_val ./ x_postdiv_s);
dx2_f = norm(test_x_postdiv_f.insert(x_input) - x_val ./ x_postdiv);
dx2_fs = norm(test_x_postdiv_fs.insert(x_input) - x_val ./ x_postdiv_s);

% add to output array
test_array = [test_array, dx2, dx2_s, dx2_f, dx2_fs];
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Array-Post-Div', 'Matrix-Scalar Array-Post-Div', ...
    'FMat-Matrix Array-Post-Div  ', 'FMat-Scalar Array-Post-Div  '};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(40, 1);
    else
        y = rome_model_var(1, 40);
    end    
    y_input = rand(size(y));
    y_postdiv = rand(size(y));
    y_postdiv_s = rand;

    % perform computation
    test_y_postdiv   = y ./ y_postdiv;
    test_y_postdiv_s = y ./ y_postdiv_s;

    y_val = y.insert(y_input);
    dy2 = norm(test_y_postdiv.insert(y_input) - y_val ./ y_postdiv);
    dy2_s = norm(test_y_postdiv_s.insert(y_input) - y_val ./ y_postdiv_s);

    % add to output array
    test_array = [test_array, dy2, dy2_s];
    descr_array = {descr_array{:}, ...
        'Vector-Vector Array-Post-Div', 'Vector-Scalar Array-Post-Div'};
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));
y_postdiv = rand(300, 30);
y_postdiv_v = rand(1, 340);
y_postdiv_s = rand;

% perform computation
test_y_postdiv   = y ./ y_postdiv;
test_y_postdiv_v = y ./ y_postdiv_v;
test_y_postdiv_s = y ./ y_postdiv_s;

y_val = y.insert(y_input);
dy2  = norm(test_y_postdiv.insert(y_input) - y_val ./ y_postdiv);
dy2_v = norm(test_y_postdiv_v.insert(y_input) - y_val ./ y_postdiv_v);
dy2_s = norm(test_y_postdiv_s.insert(y_input) - y_val ./ y_postdiv_s);

% add to output array
test_array = [test_array, dy2, dy2_v, dy2_s];
descr_array = {descr_array{:}, ...
    'Scalar-Matrix Array-Post-Div', ...
    'Scalar-Vector Array-Post-Div', ...
    'Scalar-Scalar Array-Post-Div'};


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
