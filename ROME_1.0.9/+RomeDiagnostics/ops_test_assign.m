% ops_test_assign.m
% Script to test subscripted assignment

% 10. Testing subscripted assignment
% -----------------------------------
disp ('10. Testing Subscripted Assignment...');

% For N-D Arrays
Niter = 20;
N = 8;
M = 15;
x_orig = rome_model_var(N, N, N, N);
y_orig = rome_model_var(M, M, M, M);
Ndims  = numel(size(x_orig));

% Define Inputs
x_input = rand(size(x_orig));
y_input = rand(size(y_orig));
all_input = [x_input(:); y_input(:)];

% perform subscripted assignment
for ii = 1:Niter
    x = x_orig;     % current x
    y = y_orig;     % current y    
    x_val = x_input; % test value of x
    y_val = y_input; % test value of y

    Sx = unidrnd(N, [Ndims, 1]);   % startX
    L  = unidrnd(N, [Ndims, 1]);   % length of subscripting
    Sy = unidrnd(M-L+1, [Ndims, 1]); % startY
    
    % get indices
    indX = cell(1, Ndims);
    indY = cell(1, Ndims);
    for jj = 1:Ndims
        indX{jj} = (Sx(jj)-1) + (1:L(jj));
        indY{jj} = (Sy(jj)-1) + (1:L(jj));
    end
    
    % convert some of the indices to ':'
    if(ii <= 2*Ndims)
        indX{mod(ii-1, Ndims) + 1} = ':';
        indY{mod(ii-1, Ndims) + 1} = (unidrnd(M-N)-1) + (1:N);
    end
    
    % do assignment
    x(indX{:}) = y(indY{:});    % rome_var
    x_val(indX{:}) = y_val(indY{:});    % numeric
    
    % compute error
    dx = x_val - x.insert(all_input(1:x.NumMappedVars));
    dx = norm(dx(:));
    
    % insert into results array
    test_array = [test_array, dx];    
    if(any(Sx + L > N))
        descr_msg  = sprintf('N-Dims SubsAsgn %.3d (E)', ii);   % with expansion
    else
        descr_msg  = sprintf('N-Dims SubsAsgn %.3d    ', ii);
    end    
    descr_array = {descr_array{:}, descr_msg};
end

% For Matrices
Niter = 20;
N = 20;
M = 30;
x_orig = rome_model_var(N, N);
y_orig = rome_model_var(M, M);
Ndims  = numel(size(x_orig));

% Define Inputs
x_input = rand(size(x_orig));
y_input = rand(size(y_orig));
all_input = [x_input(:); y_input(:)];

% perform subscripted assignment
for ii = 1:Niter
    x = x_orig;     % current x
    y = y_orig;     % current y
    x_val = x_input; % test value of x
    y_val = y_input; % test value of y

    Sx = unidrnd(N, [2, 1]);   % startX
    L  = unidrnd(N, [2, 1]);   % length of subscripting
    Sy = unidrnd(M-L+1, [2, 1]); % startY
    
    % get indices
    indX = cell(1, Ndims);
    indY = cell(1, Ndims);
    for jj = 1:Ndims
        indX{jj} = (Sx(jj)-1) + (1:L(jj));
        indY{jj} = (Sy(jj)-1) + (1:L(jj));
    end
    
    % convert some of the indices to ':'
    if(ii <= 2*Ndims)
        indX{mod(ii-1, Ndims) + 1} = ':';
        indY{mod(ii-1, Ndims) + 1} = (unidrnd(M-N)-1) + (1:N);
    end
    
    % do assignment
    x(indX{:}) = y(indY{:});    % rome_var
    x_val(indX{:}) = y_val(indY{:});    % numeric
    
    
    % compute error
    dx = x_val - x.insert(all_input(1:x.NumMappedVars));
    dx = norm(dx(:));
    
    % insert into results array
    test_array = [test_array, dx];    
    if(any(Sx + L > N))
        descr_msg  = sprintf('Matrix SubsAsgn %.3d (E)', ii);   % with expansion
    else
        descr_msg  = sprintf('Matrix SubsAsgn %.3d    ', ii);
    end    
    descr_array = {descr_array{:}, descr_msg};
end

% For Vectors
Niter = 10;
N = 20;
M = 30;
x_orig = rome_model_var(N, 1);
y_orig = rome_model_var(M,  1);

% Define Inputs
x_input = rand(size(x_orig));
y_input = rand(size(y_orig));
all_input = [x_input(:); y_input(:)];

% perform subscripted assignment
for ii = 1:2*Niter
    if(ii <= Niter)
        x = x_orig;     % current x   (COL vector)
        y = y_orig;     % current y
        x_val = x_input; % test value of x
        y_val = y_input; % test value of y
    else
        x = x_orig';     % current x  (ROW vector)
        y = y_orig';     % current y
        x_val = x_input'; % test value of x
        y_val = y_input'; % test value of y
    end

    Sx = unidrnd(N);   % startX
    L  = unidrnd(N);   % length of subscripting
    Sy = unidrnd(M-L+1); % startY
    
    % do subscripted assignment
    x((Sx-1) + (1:L)) = y((Sy-1) + (1:L));  % for rome_var
    x_val((Sx-1) + (1:L)) = y_val((Sy-1) + (1:L));  % for numeric
    
    % compute error
    dx = norm(x_val - x.insert(all_input(1:x.NumMappedVars)));
    
    % insert into results array
    test_array = [test_array, dx];    
    if(Sx + L > N)
        descr_msg  = sprintf('Vector SubsAsgn %.3d (E)', ii);   % with expansion
    else
        descr_msg  = sprintf('Vector SubsAsgn %.3d    ', ii);
    end
    
    descr_array = {descr_array{:}, descr_msg};
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
