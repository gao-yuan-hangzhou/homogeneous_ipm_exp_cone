% ops_test_transpose.m
% Script to test transposition. 
% Test real and complex transpose. Currently the same, but might want to
% distinguish them in the future

% 6. Testing transpose
% -----------------------------
disp ('6. Testing Transpose ...');

% when variable is a matrix
x = rome_model_var(40, 60);
x_input = rand(size(x));

% perform computation
test_x_real_tpose_1 = x.';
test_x_herm_tpose_2 = x';

% compute errors
x_val = x.insert(x_input);
dx1 = norm(test_x_real_tpose_1.insert(x_input) - x_val.');
dx2 = norm(test_x_herm_tpose_2.insert(x_input) - x_val');

% sub to output array
test_array = [test_array, dx1, dx2];
descr_array = {descr_array{:}, 'Matrix Real Transpose', 'Matrix Herm Transpose'};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(100, 1);
    else
        y = rome_model_var(1, 58);
    end    
    y_input = rand(size(y));
        
    % perform computation
    test_y_real_tpose_1    = y.';
    test_y_herm_tpose_2    = y';
    
    y_val = y.insert(y_input);
    dy1 = norm(test_y_real_tpose_1.insert(y_input) - y_val.');
    dy2 = norm(test_y_herm_tpose_2.insert(y_input) - y_val');
 
    % sub to output array
    test_array = [test_array, dy1, dy2]; 
    descr_array = {descr_array{:}, 'Vector Real Transpose', 'Vector Herm Transpose'}; 
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));

% perform computation
test_y_real_tpose_1  = y.';
test_y_herm_tpose_2  = y';

y_val = y.insert(y_input);
dy1 = norm(test_y_real_tpose_1.insert(y_input) - y_val.');
dy2 = norm(test_y_herm_tpose_2.insert(y_input) - y_val');

% sub to output array
test_array = [test_array, dy1, dy2];
descr_array = {descr_array{:}, 'Scalar Real Transpose', 'Scalar Herm Transpose'};


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
