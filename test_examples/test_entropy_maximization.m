addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;

% Whether to ensure feasibility
ensure_feasible = true;

% min sum(d(j)*x(j)*log(x(j)), j = 1:N)
% Decision variables: x = (x(1), ..., x(N))
% Constraints: Ax=b, x>=0

% The above problem can be formulated as
% min -d'*u
% Decision variables: [u(j); v(j); x(j)], j = 1, ..., N
% Constraints: Ax = b, v(j) = 1, [u(j); v(j); x(j)] in K_exp, j = 1, ..., N

disp('Constructing a random (feasible) instance...');
N = 1000; mA = 500;
d = abs(randn(N,1));
A = sprandn(mA, N, 0.05);
b = randn(mA, 1);

if ensure_feasible
    x0 = abs(randn(N,1));
    b = A*x0;
end

disp('Constructing the standard form...');
A_tilde = sparse([]); b_tilde = [];
% Code up the lines of A_tilde and b_tilde corresponding to Ax=b
for ii = 1:mA
    for j = 1:N
        % Coefficient of x(j) in constraint ii
        row = ii; 
        idx_xj = 3*j;
        A_tilde(row, idx_xj) = A(ii,j);
    end
end
b_tilde(1:mA,1) = b;

% Code up the lines of A_tilde and b_tilde correspoinding to v(j) = 1, j = 1, ..., N
for j = 1:N
    row = mA + j;
    idx_vj = 3*j-1;
    A_tilde(row, idx_vj) = 1;
    b_tilde(row) = 1;
end

% Code up vector c
c_tilde = zeros(3*N,1);
for j = 1:N
    idx_uj = 3*j-2;
    c_tilde(idx_uj) = -d(j);
end

% Assign input argument values
blk{1,1} = 'e'; blk{1,2} = 3*ones(N,1); A_cell{1} = A_tilde; c_cell{1} = c_tilde; b_input = b_tilde;

% Solve the instance
disp('Calling the solvers...');
[pd_obj, xre, yre, zre, info] = hsd_lqeu_Schur(blk, A_cell, c_cell, b_input, 1e-8);
[pd_obj, xre, yre, zre, info] = hsd_lqeu(blk, A_cell, c_cell, b_input, 1e-8);

% Retrive the optimal decision variables
x_sol = zeros(N,1);
u_sol = zeros(N,1);
v_sol = zeros(N,1);
for j = 1:N
    x_sol(j) = xre{1}(3*j);
    u_sol(j) = xre{1}(3*j-2);
    v_sol(j) = xre{1}(3*j-1);
end

obj1 = -d'*u_sol; 
obj2 = d'*(x_sol .* log(x_sol));
display(' ');
disp('Check the values of <c,xre> and sum(d(j)*x_sol(j)log(x_sol(j))):');
disp([obj1, obj2]);