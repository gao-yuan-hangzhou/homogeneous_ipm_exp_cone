% ops_test_plus.m
% Script to test addition

% 4. Testing addition
% -----------------------------
disp ('4. Testing Addition ...');

% when variable is a matrix
x = rome_model_var(50, 37);
x_input = rand(size(x));
x_preadd = rand(size(x));
x_postadd = rand(size(x));
x_preadd_s = rand;
x_postadd_s = rand;

w = rome_model_var(size(x));
w_input = rand(size(w));

% perform computation
test_x_preadd = x_preadd + x;
test_x_postadd = x + x_postadd;
test_x_preadd_s = x_preadd_s + x;
test_x_postadd_s = x + x_postadd_s;
test_x_preadd_f = x_preadd + full(x);
test_x_postadd_f = full(x) + x_postadd;
test_x_composite_1 = (x_preadd + x) + x_postadd;
test_x_composite_2 = x_preadd + (x + x_postadd);
test_xw_dbladd = x + w;
test_wx_dbladd = w + x;

% compute errors
x_val = x.insert(x_input);
dx1 = norm(test_x_preadd.insert(x_input) - (x_preadd + x_val));
dx2 = norm(test_x_postadd.insert(x_input) - (x_val + x_postadd));
dx1_s = norm(test_x_preadd_s.insert(x_input) - (x_preadd_s + x_val));
dx2_s = norm(test_x_postadd_s.insert(x_input) - (x_val + x_postadd_s));
dx1_f = norm(test_x_preadd_f.insert(x_input) - (x_preadd + x_val));
dx2_f = norm(test_x_postadd_f.insert(x_input) - (x_val + x_postadd));
dx1_c = norm(test_x_composite_1.insert(x_input) - ((x_preadd + x_val) + x_postadd));
dx2_c = norm(test_x_composite_2.insert(x_input) - (x_preadd + (x_val + x_postadd)));

w_val = w.insert(w_input);
dxw_1 = norm(test_xw_dbladd.insert([x_input(:); w_input(:)]) - (x_val + w_val));
dwx_2 = norm(test_wx_dbladd.insert([x_input(:); w_input(:)]) - (w_val + x_val));

% add to output array
test_array = [test_array, dx1, dx2, dx1_s, dx2_s, dx1_f, dx2_f, dx1_c, dx2_c, dxw_1, dwx_2];
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Pre-Add ', 'Matrix-Matrix Post-Add', ...
    'Matrix-Scalar Pre-Add ', 'Matrix-Scalar Post-Add', ...
    'FMat-Matrix Pre-Add   ', 'FMat-Matrix Post-Add  ', ...
    'Composite Array Add 1 ', 'Composite Array Add 2 ', ...
    'Matrix-Matrix Dbl-Add1', 'Matrix-Matrix Dbl-Add2'};

% when variable is a vector
for ii = 1:2
    if(ii == 0)
        y = rome_model_var(80, 1);
    else
        y = rome_model_var(1, 74);
    end    
    y_input = rand(size(y));
    y_preadd = rand(size(y));
    y_postadd = rand(size(y));
    y_preadd_s = rand;
    y_postadd_s = rand;
        
    z = rome_model_var(size(y));
    z.BiAffineMap = full(rand(size(z.BiAffineMap)));
    z_input = rand(size(z));

    % perform computation
    test_y_preadd    = y_preadd + y;
    test_y_postadd   = y + y_postadd;
    test_y_preadd_s  = y_preadd_s + y;
    test_y_postadd_s = y + y_postadd_s;
    test_yz_dbladd   = y + z;
    test_zy_dbladd   = z + y;

    y_val = y.insert(y_input);
    dy1 = norm(test_y_preadd.insert(y_input) - (y_preadd + y_val));
    dy2 = norm(test_y_postadd.insert(y_input) - (y_val + y_postadd));
    dy1_s = norm(test_y_preadd_s.insert(y_input) - (y_preadd_s + y_val));
    dy2_s = norm(test_y_postadd_s.insert(y_input) - (y_val + y_postadd_s));
    
    z_val = z.insert(z_input);
    dyz_1 = norm(test_yz_dbladd.insert([y_input(:); z_input(:)]) - (y_val + z_val));
    dzy_1 = norm(test_zy_dbladd.insert([y_input(:); z_input(:)]) - (z_val + y_val));

    % add to output array
    test_array = [test_array, dy1, dy2, dy1_s, dy2_s, dyz_1, dzy_1];
    descr_array = {descr_array{:}, ...
        'Vector-Vector Pre-Add ', 'Vector-Vector Post-Add', ...
        'Vector-Scalar Pre-Add ', 'Vector-Scalar Post-Add', ...
        'Vector-Vector Dbl-Add1', 'Vector-Vector Dbl-Add2'};
end

% when variable is a scalar
y = rome_model_var;
y_input = rand(size(y));
y_preadd = rand(700, 550);
y_postadd = rand(300, 30);
y_preadd_v = rand(800, 1);
y_postadd_v = rand(1, 340);
y_preadd_s = rand;
y_postadd_s = rand;

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
test_y_preadd    = y_preadd + y;
test_y_postadd   = y + y_postadd;
test_y_preadd_v  = y_preadd_v + y;
test_y_postadd_v = y + y_postadd_v;
test_y_preadd_s  = y_preadd_s + y;
test_y_postadd_s = y + y_postadd_s;

% double sums, scalar
test_yz_dbladd_s   = y + z_s;
test_zy_dbladd_s   = z_s + y;

% double sums, vector
test_yz_dbladd_v1   = y + z_v1;
test_zy_dbladd_v1   = z_v1 + y;
test_yz_dbladd_v2   = y + z_v2;
test_zy_dbladd_v2   = z_v2 + y;

% double sums, matrix
test_yz_dbladd_m   = y + z_m;
test_zy_dbladd_m   = z_m + y;

% compute errors
y_val = y.insert(y_input);
dy1  = norm(test_y_preadd.insert(y_input) - (y_preadd + y_val));
dy2  = norm(test_y_postadd.insert(y_input) - (y_val + y_postadd));
dy1_v = norm(test_y_preadd_v.insert(y_input) - (y_preadd_v + y_val));
dy2_v = norm(test_y_postadd_v.insert(y_input) - (y_val + y_postadd_v));
dy1_s = norm(test_y_preadd_s.insert(y_input) - (y_preadd_s + y_val));
dy2_s = norm(test_y_postadd_s.insert(y_input) - (y_val + y_postadd_s));

z_s_val = z_s.insert(z_s_input);
dyz_s = norm(test_yz_dbladd_s.insert([y_input; z_s_input]) - (y_val + z_s_val));
dzy_s = norm(test_zy_dbladd_s.insert([y_input; z_s_input]) - (z_s_val + y_val));

z_v1_val = z_v1.insert(z_v1_input);
dyz_v1 = norm(test_yz_dbladd_v1.insert([y_input; zeros1; z_v1_input]) - (y_val + z_v1_val));
dzy_v1 = norm(test_zy_dbladd_v1.insert([y_input; zeros1; z_v1_input]) - (z_v1_val + y_val));

z_v2_val = z_v2.insert(z_v2_input);
dyz_v2 = norm(test_yz_dbladd_v2.insert([y_input; zeros2; z_v2_input(:)]) - (y_val + z_v2_val));
dzy_v2 = norm(test_zy_dbladd_v2.insert([y_input; zeros2; z_v2_input(:)]) - (z_v2_val + y_val));

z_m_val  = z_m.insert(z_m_input);
dyz_m = norm(test_yz_dbladd_m.insert([y_input; zeros3; z_m_input(:)]) - (y_val + z_m_val));
dzy_m = norm(test_zy_dbladd_m.insert([y_input; zeros3; z_m_input(:)]) - (z_m_val + y_val));

% add to output array
test_array = [test_array, dy1, dy2, dy1_v, dy2_v, dy1_s, dy2_s, ...
                dyz_s, dzy_s, dyz_v1, dzy_v1, dyz_v2, dzy_v2, dyz_m, dzy_m];
descr_array = {descr_array{:}, ...
    'Scalar-Matrix Pre-Add ', 'Scalar-Matrix Post-Add', ...
    'Scalar-Vector Pre-Add ', 'Scalar-Vector Post-Add', ...
    'Scalar-Scalar Pre-Add ', 'Scalar-Scalar Post-Add', ...
    'Scalar-Matrix Dbl-Add1', 'Scalar-Matrix Dbl-Add2', ...
    'Scalar-Vector Dbl-Add1', 'Scalar-Vector Dbl-Add2', ...
    'Scalar-Vector Dbl-Add3', 'Scalar-Vector Dbl-Add4', ...
    'Scalar-Scalar Dbl-Add1', 'Scalar-Scalar Dbl-Add2'};


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
