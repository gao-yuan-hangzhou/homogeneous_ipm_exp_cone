% Reference: http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf

% The following code demonstrates the dense column handling technique for the Schur complement
% equation from the linear system for the search direction in the Homogeneous model.

m = 7; 
b = randn(m,1); 
b_bar = randn(m,1);
g_bar = randn();
B_hat = [zeros(m,m), -b, b_bar; b', 0, g_bar; -b_bar', -g_bar, 0];
h_hat = randn(m+2,1);

U1 = [-b; 0; -g_bar/2];
U2 = [b_bar; g_bar/2; 0];

V1 = [zeros(m,1); 1; 0];
V2 = [zeros(m+1,1); 1];

U = [U1, U2, V1, V2];
D = [zeros(2), eye(2); -eye(2), zeros(2)];

disp(['||U*D*U'' - Bhat|| = ' num2str(norm(U*D*U'-B_hat))]);

Msp = randn(m+2, m+2);
Msp = Msp'*Msp;

M = Msp + U*D*U';

Dinv = [zeros(2), -eye(2); eye(2), zeros(2)];
E = -Dinv;

% Compute the actual solution of the equation M*yHat = h
dy_hat_actual = M\h_hat;

% Now compute the solution through the dense column handling technique
% Define lHat = D*U'*yHat and it can be shown that the system M*yHat = hHat is equivalent to 
% Mbig*[dyHat; lHat] = hBig, where
% M_big = [Msp, U; U' , E]
% h_big = [h_hat; zeros(4,1)]
% dy_big = M_big\h_big

% Now solve for dy_big using Sherman-Morrison formula
[L,p] = chol(Msp);
Msp_inv = ;
G = eye()

% dy_hat_dch = dy_big(1:m+2);

norm(dy_hat_actual - dy_hat_dch)