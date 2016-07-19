% ops_test_sum_rand.m
% Script to test summation and subscripting, using non-structured indices

% 8. Testing summation and subscripting
% ------------------------------------
disp ('8. Testing Summation & Subscripting (Non-Structured)...');

% when variable is a multilinear form
x = rome_model_var(10, 10, 5, 6, 4);
x_input = rand(size(x));
x_val = x.insert(x_input);
sz_x = size(x);
N_iter = (length(size(x))+10);
summing_dim = unidrnd(length(size(x)), N_iter, 1);

for ii = 1:N_iter
    % perform computation, sum over a random index set
%     ind_limit = unidrnd(sz_x);
    ind_limit = sz_x;
    ind = cell(1, length(sz_x));
    for jj = 1:length(sz_x)
        ind{jj} = unidrnd(sz_x(jj), [ind_limit(jj), 1]);
    end
    
    if(ii <= length(sz_x))
        % change one of the indices to ':'
        ind{ii} = ':';
    end
    
    % perform computation
    w = x(ind{:});
    if(mod(ii, 2) == 0)        
        test_sum = sum(w);
%         test_sum = sum(x(ind{:}));
        err = test_sum.insert(x_input) - sum(x_val(ind{:}));
        dx = norm(err(:));
    else
        test_sum = sum(w, summing_dim(ii));
%         test_sum = sum(x(ind{:}), summing_dim(ii));
        err = test_sum.insert(x_input) - sum(x_val(ind{:}), summing_dim(ii));
        dx = norm(err(:));
    end
    
    % compute errors
    
    test_array = [test_array, dx];
    descr_array = {descr_array{:}, sprintf('MultiLin Sum %0.3d', ii)};
end

% when variable is a matrix
x = rome_model_var(60, 75);
x_input = rand(size(x));
x_val = x.insert(x_input);
sz_x = size(x);
N_iter = length(sz_x) + 15;
summing_dim = unidrnd(length(sz_x), N_iter, 1);

for ii = 1:N_iter
    % perform computation, sum over a random index set
    ind_limit = unidrnd(sz_x);
    ind = cell(1, length(sz_x));
    for jj = 1:length(sz_x)
        ind{jj} = unidrnd(sz_x(jj), [ind_limit(jj), 1]);
    end
    
    if(ii <= length(sz_x))
        % change one of the indices to ':'
        ind{ii} = ':';
    end
    
    % perform computation
    if(mod(ii, 2) == 0)
        test_sum = sum(x(ind{:}));
        dx = norm(test_sum.insert(x_input) - sum(x_val(ind{:})));
    else
        test_sum = sum(x(ind{:}), summing_dim(ii));
        dx = norm(test_sum.insert(x_input) - sum(x_val(ind{:}), summing_dim(ii)));
    end
    
    % compute errors
    test_array = [test_array, dx];
    descr_array = {descr_array{:}, sprintf('Matrix Sum %0.3d', ii)};
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
    
    for ii = 1:20
        % perform computation, sum over a random index set
        ind_limit = unidrnd(length(y));
        ind = unidrnd(length(y), [ind_limit, 1]);

        % perform computation
        test_sum = sum(y(ind));

        % compute errors
        dy = norm(test_sum.insert(y_input) - y_val(ind));
        test_array = [test_array, dy];
        descr_array = {descr_array{:}, sprintf('Vector Sum %0.3d', ii)};
    end    
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
