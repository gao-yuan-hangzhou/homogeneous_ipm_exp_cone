% ops_test_minus.m
% Script to test subtraction

% 5. Testing subtraction
% -----------------------------
disp ('5. Testing Subtraction ...');

% when variable is a matrix
x = rome_model_var(53, 39);
x_input = rand(size(x));
x_presub = rand(size(x));
x_postsub = rand(size(x));
x_presub_s = rand;
x_postsub_s = rand;

w = rome_model_var(size(x));
w_input = rand(size(w));

% perform computation
test_x_presub = x_presub - x;
test_x_postsub = x - x_postsub;
test_x_presub_s = x_presub_s - x;
test_x_postsub_s = x - x_postsub_s;
test_x_presub_f = x_presub - full(x);
test_x_postsub_f = full(x) - x_postsub;
test_x_composite_1 = (x_presub - x) - x_postsub;
test_x_composite_2 = x_presub - (x - x_postsub);
test_xw_dblsub = x - w;
test_wx_dblsub = w - x;

% compute errors
x_val = x.insert(x_input);
dx1 = norm(test_x_presub.insert(x_input) - (x_presub - x_val));
dx2 = norm(test_x_postsub.insert(x_input) - (x_val - x_postsub));
dx1_s = norm(test_x_presub_s.insert(x_input) - (x_presub_s - x_val));
dx2_s = norm(test_x_postsub_s.insert(x_input) - (x_val - x_postsub_s));
dx1_f = norm(test_x_presub_f.insert(x_input) - (x_presub - x_val));
dx2_f = norm(test_x_postsub_f.insert(x_input) - (x_val - x_postsub));
dx1_c = norm(test_x_composite_1.insert(x_input) - ((x_presub - x_val) - x_postsub));
dx2_c = norm(test_x_composite_2.insert(x_input) - (x_presub - (x_val - x_postsub)));

w_val = w.insert(w_input);
dxw_1 = norm(test_xw_dblsub.insert([x_input(:); w_input(:)]) - (x_val - w_val));
dwx_2 = norm(test_wx_dblsub.insert([x_input(:); w_input(:)]) - (w_val - x_val));

% sub to output array
test_array = [test_array, dx1, dx2, dx1_s, dx2_s, dx1_f, dx2_f, dx1_c, dx2_c, dxw_1, dwx_2];
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Pre-Sub ', 'Matrix-Matrix Post-Sub', ...
    'Matrix-Scalar Pre-Sub ', 'Matrix-Scalar Post-Sub', ...
    'FMat-Matrix Pre-Sub   ', 'FMat-Matrix Post-Sub  ', ...
    'Composite Array Sub 1 ', 'Composite Array Sub 2 ', ...
    'Matrix-Matrix Dbl-Sub1', 'Matrix-Matrix Dbl-Sub2'};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(80, 1);
    else
        y = rome_model_var(1, 74);
    end    
    y_input = rand(size(y));
    y_presub = rand(size(y));
    y_postsub = rand(size(y));
    y_presub_s = rand;
    y_postsub_s = rand;
        
    z = rome_model_var(size(y));
    z.BiAffineMap = full(rand(size(z.BiAffineMap)));
    z_input = rand(size(z));

    % perform computation
    test_y_presub    = y_presub - y;
    test_y_postsub   = y - y_postsub;
    test_y_presub_s  = y_presub_s - y;
    test_y_postsub_s = y - y_postsub_s;
    test_yz_dblsub   = y - z;
    test_zy_dblsub   = z - y;

    y_val = y.insert(y_input);
    dy1 = norm(test_y_presub.insert(y_input) - (y_presub - y_val));
    dy2 = norm(test_y_postsub.insert(y_input) - (y_val - y_postsub));
    dy1_s = norm(test_y_presub_s.insert(y_input) - (y_presub_s - y_val));
    dy2_s = norm(test_y_postsub_s.insert(y_input) - (y_val - y_postsub_s));
    
    z_val = z.insert(z_input);
    dyz_1 = norm(test_yz_dblsub.insert([y_input(:); z_input(:)]) - (y_val - z_val));
    dzy_1 = norm(test_zy_dblsub.insert([y_input(:); z_input(:)]) - (z_val - y_val));

    % sub to output array
    test_array = [test_array, dy1, dy2, dy1_s, dy2_s, dyz_1, dzy_1];
    descr_array = {descr_array{:}, ...
        'Vector-Vector Pre-Sub ', 'Vector-Vector Post-Sub', ...
        'Vector-Scalar Pre-Sub ', 'Vector-Scalar Post-Sub', ...
        'Vector-Vector Dbl-Sub1', 'Vector-Vector Dbl-Sub2'};
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));
y_presub = rand(700, 550);
y_postsub = rand(300, 30);
y_presub_v = rand(800, 1);
y_postsub_v = rand(1, 340);
y_presub_s = rand;
y_postsub_s = rand;

z_s  = rome_model_var;
z_v1 = rome_model_var(45, 1);
z_v2 = rome_model_var(1, 60);
z_m  = rome_model_var(80, 75);

% create zeros for insertion later
zeros1 = zeros(prod(size(z_s)), 1);
zeros2 = [zeros1; zeros(prod(size(z_v1)), 1)];
zeros3 = [zeros2; zeros(prod(size(z_v2)), 1)];
zeros4 = [zeros3; zeros(prod(size(z_m)), 1)];


z_s_input = rand(size(z_s));
z_v1_input = rand(size(z_v1));
z_v2_input = rand(size(z_v2));
z_m_input  = rand(size(z_m));

% perform computation
test_y_presub    = y_presub - y;
test_y_postsub   = y - y_postsub;
test_y_presub_v  = y_presub_v - y;
test_y_postsub_v = y - y_postsub_v;
test_y_presub_s  = y_presub_s - y;
test_y_postsub_s = y - y_postsub_s;

% double sums, scalar
test_yz_dblsub_s   = y - z_s;
test_zy_dblsub_s   = z_s - y;

% double sums, vector
test_yz_dblsub_v1   = y - z_v1;
test_zy_dblsub_v1   = z_v1 - y;
test_yz_dblsub_v2   = y - z_v2;
test_zy_dblsub_v2   = z_v2 - y;

% double sums, matrix
test_yz_dblsub_m   = y - z_m;
test_zy_dblsub_m   = z_m - y;

% compute errors
y_val = y.insert(y_input);
dy1  = norm(test_y_presub.insert(y_input) - (y_presub - y_val));
dy2  = norm(test_y_postsub.insert(y_input) - (y_val - y_postsub));
dy1_v = norm(test_y_presub_v.insert(y_input) - (y_presub_v - y_val));
dy2_v = norm(test_y_postsub_v.insert(y_input) - (y_val - y_postsub_v));
dy1_s = norm(test_y_presub_s.insert(y_input) - (y_presub_s - y_val));
dy2_s = norm(test_y_postsub_s.insert(y_input) - (y_val - y_postsub_s));

z_s_val = z_s.insert(z_s_input);
dyz_s = norm(test_yz_dblsub_s.insert([y_input; z_s_input]) - (y_val - z_s_val));
dzy_s = norm(test_zy_dblsub_s.insert([y_input; z_s_input]) - (z_s_val - y_val));

z_v1_val = z_v1.insert(z_v1_input);
dyz_v1 = norm(test_yz_dblsub_v1.insert([y_input; zeros1; z_v1_input]) - (y_val - z_v1_val));
dzy_v1 = norm(test_zy_dblsub_v1.insert([y_input; zeros1; z_v1_input]) - (z_v1_val - y_val));

z_v2_val = z_v2.insert(z_v2_input);
dyz_v2 = norm(test_yz_dblsub_v2.insert([y_input; zeros2; z_v2_input(:)]) - (y_val - z_v2_val));
dzy_v2 = norm(test_zy_dblsub_v2.insert([y_input; zeros2; z_v2_input(:)]) - (z_v2_val - y_val));

z_m_val  = z_m.insert(z_m_input);
dyz_m = norm(test_yz_dblsub_m.insert([y_input; zeros3; z_m_input(:)]) - (y_val - z_m_val));
dzy_m = norm(test_zy_dblsub_m.insert([y_input; zeros3; z_m_input(:)]) - (z_m_val - y_val));

% sub to output array
test_array = [test_array, dy1, dy2, dy1_v, dy2_v, dy1_s, dy2_s, ...
                dyz_s, dzy_s, dyz_v1, dzy_v1, dyz_v2, dzy_v2, dyz_m, dzy_m];
descr_array = {descr_array{:}, ...
    'Scalar-Matrix Pre-Sub ', 'Scalar-Matrix Post-Sub', ...
    'Scalar-Vector Pre-Sub ', 'Scalar-Vector Post-Sub', ...
    'Scalar-Scalar Pre-Sub ', 'Scalar-Scalar Post-Sub', ...
    'Scalar-Matrix Dbl-Sub1', 'Scalar-Matrix Dbl-Sub2', ...
    'Scalar-Vector Dbl-Sub1', 'Scalar-Vector Dbl-Sub2', ...
    'Scalar-Vector Dbl-Sub3', 'Scalar-Vector Dbl-Sub4', ...
    'Scalar-Scalar Dbl-Sub1', 'Scalar-Scalar Dbl-Sub2'};


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
