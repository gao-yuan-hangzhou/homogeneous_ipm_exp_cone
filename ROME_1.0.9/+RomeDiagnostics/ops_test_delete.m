% ops_test_delete.m
% Script to test subscripted deletion

% 11. Testing subscripted deletion
% -----------------------------------
disp ('11. Testing Subscripted Deletion...');

% For N-D Arrays
Niter = 20;
N = 8;
x_orig = rome_model_var(N, N, N, N);
Ndims  = numel(size(x_orig));

% Define Inputs
x_input = rand(size(x_orig));

% perform subscripted assignment
ind = cell(1, Ndims);
for ii = 1:Niter
    x     = x_orig;     % current x
    x_val = x_input; % test value of x

    L  = unidrnd(N);    % length of subscripting
    Sx = unidrnd(N-L+1);  % startX
    
    % make deletion indices
    [ind{:}] = deal(':');
    ind{mod(ii-1, Ndims) + 1} = (Sx-1) + (1:L);
    
    % do deletion
    x(ind{:}) = [];        % rome_var
    x_val(ind{:}) = [];    % numeric
    
    % compute error
    if(isempty(x_val))
        dx = 0;        
        descr_msg = sprintf('N-Dims Delete %.3d (e) ', ii);
    else
        dx = x_val - x.insert(x_input);
        dx = norm(dx(:));        
        descr_msg = sprintf('N-Dims Delete %.3d     ', ii);
    end
    
    % insert into results array
    test_array = [test_array, dx];
    
    % Message
    descr_array = {descr_array{:}, descr_msg};
end

% For Matrices
Niter = 20;
N = 8;
x_orig = rome_model_var(N, N);
Ndims  = numel(size(x_orig));

% Define Inputs
x_input = rand(size(x_orig));

% perform subscripted assignment
ind = cell(1, Ndims);
for ii = 1:Niter
    x     = x_orig;     % current x
    x_val = x_input; % test value of x

    L  = unidrnd(N);    % length of subscripting
    Sx = unidrnd(N-L+1);  % startX
    
    % make deletion indices
    [ind{:}] = deal(':');
    ind{mod(ii-1, Ndims) + 1} = (Sx-1) + (1:L);
    
    % do deletion
    x(ind{:}) = [];        % rome_var
    x_val(ind{:}) = [];    % numeric
    
     % compute error
    if(isempty(x_val))
        dx = 0;        
        descr_msg = sprintf('Matrix Delete %.3d (e) ', ii);
    else
        dx = x_val - x.insert(x_input);
        dx = norm(dx(:));        
        descr_msg = sprintf('Matrix Delete %.3d     ', ii);
    end
    
    % insert into results array
    test_array = [test_array, dx];
    
    % Message
    descr_array = {descr_array{:}, descr_msg};
end

% For Vectors
Niter = 10;
N = 8;
x_orig = rome_model_var(N, 1);

% Define Inputs
x_input = rand(size(x_orig));

% perform subscripted assignment
for ii = 1:2*Niter
    if(ii <= Niter)
        x     = x_orig;     % current x (COL Vector)
        x_val = x_input;    % test value of x
    else
        x     = x_orig';     % current x (ROW Vector)
        x_val = x_input';    % test value of x
    end

    L  = unidrnd(N);    % length of subscripting
    Sx = unidrnd(N-L+1);  % startX
    
    % make deletion indices
    ind = (Sx-1) + (1:L);
    
    % do deletion
    x(ind) = [];        % rome_var
    x_val(ind) = [];    % numeric
    
     % compute error
    if(isempty(x_val))
        dx = 0;        
        descr_msg = sprintf('Vector Delete %.3d (e) ', ii);
    else
        dx = x_val - x.insert(x_input);
        dx = norm(dx(:));        
        descr_msg = sprintf('Vector Delete %.3d     ', ii);
    end
    
    % insert into results array
    test_array = [test_array, dx];
    
    % Message
    descr_array = {descr_array{:}, descr_msg};
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
