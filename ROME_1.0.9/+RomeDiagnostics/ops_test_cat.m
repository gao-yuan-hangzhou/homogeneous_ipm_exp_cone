% ops_test_cat
% Script To Test Concatenation

% 9. Testing array multiplication
% --------------------------------
disp ('9. Testing Concatenation ...');

% Test Cat 1-3
% ----------------
% define model variables
x = rome_model_var(13, 5, 15, 3); x_input = rand(size(x));

y1 = rome_model_var(7, 5, 15, 3); y1_input = rand(size(y1));
z1 = rome_model_var(2, 5, 15, 3); z1_input = rand(size(z1));

y2 = rome_model_var(13, 7, 15, 3); y2_input = rand(size(y2));
z2 = rome_model_var(13, 2, 15, 3); z2_input = rand(size(z2));

y3 = rome_model_var(13, 5, 7, 3); y3_input = rand(size(y3));
z3 = rome_model_var(13, 5, 2, 3); z3_input = rand(size(z3));

y4 = rome_model_var(13, 5, 15, 7); y4_input = rand(size(y4));
z4 = rome_model_var(13, 5, 15, 2); z4_input = rand(size(z4));

% perform concatenation
test_cat1a = cat(1, x, y1, z1);
test_cat1b = cat(1, z1, x, y1);
test_cat1c = cat(1, y1, z1, x);
test_cat2a = cat(2, x, y2, z2);
test_cat2b = cat(2, z2, x, y2);
test_cat2c = cat(2, y2, z2, x);
test_cat3a = cat(3, x, y3, z3);
test_cat3b = cat(3, z3, x, y3);
test_cat3c = cat(3, y3, z3, x);
test_cat4a = cat(4, x, y4, z4);
test_cat4b = cat(4, z4, x, y4);
test_cat4c = cat(4, y4, z4, x);

v1_input = [x_input(:); y1_input(:); z1_input(:)];
v2_input = [x_input(:); zeros(numel(y1_input) + numel(z1_input), 1); y2_input(:); z2_input(:)];
v3_input = [x_input(:); zeros(numel(y1_input) + numel(y2_input) ...
                            + numel(z1_input) + numel(z2_input), 1); ...
                              y3_input(:); z3_input(:)];
v4_input = [x_input(:); zeros(numel(y1_input) + numel(y2_input) + numel(y3_input) ...
                            + numel(z1_input) + numel(z2_input) + numel(z3_input), 1); ...
                              y4_input(:); z4_input(:)];

% define ground truth
v1_val_a = cat(1, x_input, y1_input, z1_input);
v1_val_b = cat(1, z1_input, x_input, y1_input);
v1_val_c = cat(1, y1_input, z1_input, x_input);
v2_val_a = cat(2, x_input, y2_input, z2_input);
v2_val_b = cat(2, z2_input, x_input, y2_input);
v2_val_c = cat(2, y2_input, z2_input, x_input);
v3_val_a = cat(3, x_input, y3_input, z3_input);
v3_val_b = cat(3, z3_input, x_input, y3_input);
v3_val_c = cat(3, y3_input, z3_input, x_input);
v4_val_a = cat(4, x_input, y4_input, z4_input);
v4_val_b = cat(4, z4_input, x_input, y4_input);
v4_val_c = cat(4, y4_input, z4_input, x_input);

% find errors 
dv1a = test_cat1a.insert(v1_input) - v1_val_a;
dv1b = test_cat1b.insert(v1_input) - v1_val_b;
dv1c = test_cat1c.insert(v1_input) - v1_val_c;

dv2a = test_cat2a.insert(v2_input) - v2_val_a;
dv2b = test_cat2b.insert(v2_input) - v2_val_b;
dv2c = test_cat2c.insert(v2_input) - v2_val_c;

dv3a = test_cat3a.insert(v3_input) - v3_val_a;
dv3b = test_cat3b.insert(v3_input) - v3_val_b;
dv3c = test_cat3c.insert(v3_input) - v3_val_c;

dv4a = test_cat4a.insert(v4_input) - v4_val_a;
dv4b = test_cat4b.insert(v4_input) - v4_val_b;
dv4c = test_cat4c.insert(v4_input) - v4_val_c;

% compute norms
dv1a = norm(dv1a(:));
dv1b = norm(dv1b(:));
dv1c = norm(dv1c(:));
dv2a = norm(dv2a(:));
dv2b = norm(dv2b(:));
dv2c = norm(dv2c(:));
dv3a = norm(dv3a(:));
dv3b = norm(dv3b(:));
dv3c = norm(dv3c(:));
dv4a = norm(dv4a(:));
dv4b = norm(dv4b(:));
dv4c = norm(dv4c(:));

% add to output array
test_array = [test_array, dv1a, dv1b, dv1c, dv2a, dv2b, dv2c, ...
                          dv3a, dv3b, dv3c, dv4a, dv4b, dv4c]; 
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Cat 1A', 'Matrix-Matrix Cat 1B', 'Matrix-Matrix Cat 1C', ...
    'Matrix-Matrix Cat 2A', 'Matrix-Matrix Cat 2B', 'Matrix-Matrix Cat 2C', ...
    'Matrix-Matrix Cat 3A', 'Matrix-Matrix Cat 3B', 'Matrix-Matrix Cat 3C', ...
    'Matrix-Matrix Cat 4A', 'Matrix-Matrix Cat 4B', 'Matrix-Matrix Cat 4C', ...
    };





% Test VertCat
% -------------
x = rome_model_var(37, 50);
x_input = rand(size(x));

y = rome_model_var(123, 50);
y_input = rand(size(y));

test_vertcat_1 = [x; y];
test_vertcat_2 = [y; x];

v_input = [x_input(:); y_input(:)];
v1_val = [x_input; y_input];
v2_val = [y_input; x_input];

dv1 = norm(test_vertcat_1.insert(v_input) - v1_val, 'fro');
dv2 = norm(test_vertcat_2.insert(v_input) - v2_val, 'fro');

% add to output array
test_array = [test_array, dv1, dv2]; 
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Vertcat 1', 'Matrix-Matrix Vertcat 2', ...
    };

% Test HorzCat
% ------------
x = rome_model_var(60, 85);
x_input = rand(size(x));

y = rome_model_var(60, 36);
y_input = rand(size(y));

test_horzcat_1 = [x, y];
test_horzcat_2 = [y, x];

h_input = [x_input(:); y_input(:)];
h1_val = [x_input, y_input];
h2_val = [y_input, x_input];

dh1 = norm(test_horzcat_1.insert(h_input) - h1_val, 'fro');
dh2 = norm(test_horzcat_2.insert(h_input) - h2_val, 'fro');

% add to output array
test_array = [test_array, dh1, dh2]; 
descr_array = {descr_array{:}, ...
    'Matrix-Matrix Horzcat 1', 'Matrix-Matrix Horzcat 2', ...
    };

% Test Reshape
% ------------
x = rome_model_var(60, 30);
x_input = rand(size(x));

tot_size_x = x.TotalSize;
sz_factor = factor(tot_size_x);
dy_arr = zeros(1, numel(sz_factor));
for ii = 1:numel(sz_factor)
    new_sz = [sz_factor(ii), floor(tot_size_x ./ sz_factor(ii))];
    y = reshape(x, new_sz);
    dy_arr(ii) = norm(y.insert(x_input(:)) - reshape(x_input, new_sz), 'fro');
end

% add to output array
test_array = [test_array, dy_arr];

for ii = 1:numel(sz_factor)
    descr_array = {descr_array{:}, sprintf('Matrix-Matrix Reshape %d', ii)};
end

% Test repmat
% --------------
x = rome_model_var(40, 50);
x_input = rand(size(x));
pattern_array = {[3, 4]; ...
                 [4, 1]; ...
                 [1, 2]; ...
                 [4, 2, 3]; ...
                 [2, 3, 1]; ...
                 [1, 2, 3]; ...
                };
dy_arr = zeros(1, size(pattern_array, 1));
for ii = 1:size(pattern_array, 1)
    y = repmat(x, pattern_array{ii});
    dy =y.insert(x_input(:)) - repmat(x_input, pattern_array{ii});
    dy_arr(ii) = norm(dy(:), 'fro');
end

% add to output array
test_array = [test_array, dy_arr];

for ii = 1:size(pattern_array, 1)
    descr_array = {descr_array{:}, sprintf('Matrix-Matrix Repmat %d', ii)};
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
