% Reference: http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf

% The following code demonstrates the dense column handling technique for the Schur complement
% equation from the linear system for the search direction in the Homogeneous model.

m = 3000; 
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

disp(['||U*D*U'' - Bhat|| = ' num2str(norm(U*D*U'-B_hat))]);

Msp = sprandn(m+2, m+2, 1/m);
Msp = Msp'*Msp + spdiags(abs(randn(m+2)), 0, m+2, m+2);
% keyboard;

M = Msp + U*D*U';

D_inv = [zeros(2), -eye(2); eye(2), zeros(2)];

% Compute the actual solution of the equation M*yHat = h
dy_hat_actual = M\h_hat;

% Now compute the solution through the dense column handling technique
% Define lHat = D*U'*yHat and it can be shown that the system M*yHat = hHat is equivalent to 
% Mbig*[dyHat; lHat] = hBig, where
% M_big = [Msp, U; U' , E]
% h_big = [h_hat; zeros(4,1)]
% dy_big = M_big\h_big

% Now solve for dy_hat using Sherman-Morrison formula
[R,p] = chol(Msp, 'lower'); R = R';
% Calculate inv(Msp)*h_hat
w_hat = R\(R'\h_hat);

% Calculate the 4-by-4 G = D_inv + U'*inv(Msp)*U and then inv(G)
M_inv_U = R\(R'\U);
G = U'*M_inv_U;
G = D_inv + G;

% Calculate P = eye(m+2) - inv(Msp)*U*G_inv*U' and then y_hat = P*w_hat
P = G\U'; 
P = M_inv_U*P;
P = eye(m+2) - P;
dy_hat = P*w_hat;

% dy_hat_dch = dy_big(1:m+2);

norm(dy_hat_actual - dy_hat)/norm(dy_hat)