% Reference: http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf

% The following code demonstrates the dense column handling technique for the Schur complement
% equation from the linear system for the search direction in the Homogeneous model.

m = 5000;
b = randn(m,1); 
b_bar = randn(m,1);
g_bar = randn();
B_hat = [sparse(m,m), -b, b_bar; b', 0, g_bar; -b_bar', -g_bar, 0];
h_hat = randn(m+2,1);

U1 = [-b; 0; -g_bar/2];
U2 = [b_bar; g_bar/2; 0];

V1 = [zeros(m,1); 1; 0];
V2 = [zeros(m+1,1); 1];

U = [U1, U2, V1, V2];
D = [zeros(2), eye(2); -eye(2), zeros(2)];

% disp(['||U*D*U'' - Bhat|| = ' num2str(norm(U*D*U'-B_hat), Inf)]);

disp('Generating random positive definite Msp...');
Msp = sprandn(m+2, m+2, 9/m);
Msp = Msp'*Msp + spdiags(abs(randn(m+2)), 0, m+2, m+2);
disp('Done...');
M = Msp + U*D*U';
% keyboard;
D_inv = [zeros(2), -eye(2); eye(2), zeros(2)];

% Compute the actual solution of the equation M*yHat = h
tic;
dy_hat_actual = M\h_hat;
toc;

% Now compute the solution through the dense column handling technique
% Define lHat = D*U'*yHat and it can be shown that the system M*yHat = hHat is equivalent to 
% Mbig*[dyHat; lHat] = hBig, where
% M_big = [Msp, U; U' , E]
% h_big = [h_hat; zeros(4,1)]
% dy_big = M_big\h_big

tic;
% Now solve for dy_hat using Sherman-Morrison formula
% [R,p] = chol(Msp);
% Calculate inv(Msp)*h_hat
Msp_inv_U_and_h_hat = Msp\[U, h_hat]; % Matlab will use sparse Cholesky factorization for this step
Msp_inv_U = Msp_inv_U_and_h_hat(:, 1:end-1);
Msp_inv_h_hat = Msp_inv_U_and_h_hat(:, end);
%Msp_inv_h_hat = Msp\h_hat;

% Calculate the 4-by-4 G = D_inv + U'*inv(Msp)*U and then 
% Calculate inv(G)
G = D_inv + U'*Msp_inv_U;

% Implicitly calculate P = eye(m+2) - inv(Msp)*U*G_inv*U' and then y_hat = P*w_hat
% P = eye(m+2) - M_inv_U*(G\U');
dy_hat = Msp_inv_h_hat - Msp_inv_U*(G\(U'*Msp_inv_h_hat));
% dy_hat_dch = dy_big(1:m+2);
toc;

lu(G);

tic;
L = chol(Msp);
Msp_inv_chol = @(X) L\(L'\X);
Msp_inv_U = Msp_inv_chol(U);
Msp_inv_h_hat = Msp_inv_chol(h_hat);
G = D_inv + U'*Msp_inv_U;
[UG, LG] = lu(G);
precond_func = @(rhs) precond_M_inv(rhs, U, L, LG, UG, Msp_inv_U);

% Compute precond(M) and precond(h_hat)
dy_hat_bicgstab = bicgstab(M, h_hat, [], [], precond_func);
toc;

disp('Comparison of the solutions y_hat_lu, y_hat_chol, y_hat_bicgs:')