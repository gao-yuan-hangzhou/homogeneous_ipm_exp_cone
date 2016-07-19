% ops_test_sum.m
% Script to test summation and subscripting

% 7. Testing summation and subscripting
% ------------------------------------
disp ('7. Testing Summation & Subscripting ...');

% when variable is a multilinear form
x = rome_model_var(10, 10, 5, 6, 4, 4);
x_input = rand(size(x));
x_val = x.insert(x_input);
sz_x = size(x);
rep_cycle = max(sz_x) - 1;
N_iter = (length(size(x)) + rep_cycle);
summing_dim = unidrnd(length(size(x)), N_iter, 1);

for ii = 1:N_iter
    % perform computation, sum over a random index set
    ind = cell(1, length(sz_x));
    cur_rep = mod(ii, rep_cycle) + 1;
    
    % sum over indices
    for jj = 1:length(sz_x)        
        ind{jj} = (1:cur_rep:sz_x(jj))';
    end
    
    if(ii <= length(sz_x))
        % change one of the indices to ':'
        ind{ii} = ':';
    end
    
    % perform computation
    x_subset = x(ind{:});
    err = x_subset.insert(x_input) - x_val(ind{:});
    
    if(mod(ii, 2) == 0)
        test_sum = sum(x_subset);
        test_diff = diff(x_subset);
        
        sum_err = test_sum.insert(x_input) - sum(x_val(ind{:}));
        diff_err = test_diff.insert(x_input) - diff(x_val(ind{:}));        
    else
        test_sum = sum(x_subset, summing_dim(ii));
        test_diff = diff(x_subset, [], summing_dim(ii));
        
        sum_err = test_sum.insert(x_input) - sum(x_val(ind{:}), summing_dim(ii));
        
        if(~isempty(test_diff))
            diff_err = test_diff.insert(x_input) - diff(x_val(ind{:}), [], summing_dim(ii));
        else
            diff_err = 0 - diff(x_val(ind{:}), [], summing_dim(ii));
        end      
    end
    
    % Compute actual errors
    dx = norm(err(:));    
    dsx = norm(sum_err(:));
    ddx = norm(diff_err(:));
    
    % put into results array
    test_array = [test_array, dx, dsx, ddx];
    descr_array = {descr_array{:}, ...
        sprintf('MultiLin Subsref %0.3d', ii), ...
        sprintf('MultiLin Summing %0.3d', ii), ...
        sprintf('MultiLin Diffing %0.3d', ii), ...
        };
end

% when variable is a matrix
x = rome_model_var(60, 75);
x_input = rand(size(x));
x_val = x.insert(x_input);
sz_x = size(x);
rep_cycle = max(sz_x) - 1;
N_iter = (length(size(x)) + rep_cycle);
summing_dim = unidrnd(length(sz_x), N_iter, 1);

for ii = 1:N_iter
    % perform computation, sum over a random index set
    ind = cell(1, length(sz_x));
    for jj = 1:length(sz_x)
        ind{jj} = (1:cur_rep:sz_x(jj))';
    end
    
    if(ii <= length(sz_x))
        % change one of the indices to ':'
        ind{ii} = ':';
    end
    
    % perform computation
    x_subset = x(ind{:});
    err = x_subset.insert(x_input) - x_val(ind{:});
    dx = norm(err(:));
    
    if(mod(ii, 2) == 0)
        test_sum = sum(x_subset);
        test_diff = diff(x_subset);
        
        dsx = norm(test_sum.insert(x_input) - sum(x_val(ind{:})));
        ddx = norm(test_diff.insert(x_input) - diff(x_val(ind{:})));
    else
        test_sum = sum(x_subset, summing_dim(ii));
        test_diff = diff(x_subset, [], summing_dim(ii));
        
        dsx = norm(test_sum.insert(x_input) - sum(x_val(ind{:}), summing_dim(ii)));
        if(~isempty(test_diff))
            ddx = norm(test_diff.insert(x_input) - diff(x_val(ind{:}), [], summing_dim(ii)));
        else
            ddx = norm(diff(x_val(ind{:}), [], summing_dim(ii)));
        end
    end
    
    % compute errors
    test_array = [test_array, dx, dsx, ddx];
    descr_array = {descr_array{:}, ...
        sprintf('Matrix Subsref %0.3d', ii), ...
        sprintf('Matrix Summing %0.3d', ii), ...
        sprintf('Matrix Diffing %0.3d', ii), ...
        };
end

% when variable is a vector
for kk = 1:2
    if(kk == 0)
        y = rome_model_var(80, 1);
    else
        y = rome_model_var(1, 74);
    end
    
    % make input
    y_input = rand(size(y));
    y_val = y.insert(y_input);
    
    for ii = 1:50
        % mirrors the input
        y_input_mirror = y_input;
        
        % perform computation, sum over a random index set
        ind_limit = unidrnd(length(y));
        ind = unidrnd(length(y), [ind_limit, 1]);

        % compute error in subscripting
        y_subset = y(ind);
        
        % work-around for contiguous indexing
        if(all(diff(ind) == 1))
            y_input_mirror = y_input(ind);
        end
        err = y_subset.insert(y_input_mirror) - y_val(ind);
        dy = norm(err(:));
        
        % perform computation
        test_sum = sum(y_subset);
        test_diff = diff(y_subset);

        % compute errors
        dsy = norm(test_sum.insert(y_input_mirror) - sum(y_val(ind)));
        if(~isempty(test_diff))
            ddy = norm(test_diff.insert(y_input_mirror) - diff(y_val(ind)));
        else
            ddy = norm(diff(y_val(ind)));
        end        
        
        test_array = [test_array, dy, dsy, ddy];
        descr_array = {descr_array{:}, ...
            sprintf('Vector Subsref %0.3d', ii), ... 
            sprintf('Vector Summing %0.3d', ii), ...
            sprintf('Vector Diffing %0.3d', ii), ...
            };
    end    
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
