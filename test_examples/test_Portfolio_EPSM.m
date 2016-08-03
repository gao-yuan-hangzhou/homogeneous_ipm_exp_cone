addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;
disp('Portfolio Allocation under the Entropic Satisficing Measure');
% Portfolio_EPSM_example
% See http://robustopt.com/examples/portfolio_epsm/Portfolio%20EPSM.pdf for further details

% The original maximization problem can be transformed into a convex minimization problem:
% min theta 
% s.t. theta = p(i)s(i)+q(i)t(i)
%      theta >= 0
%      w(1)+...+w(n) = 1, w >= 0, i = 1,...,n
%      z(1)+...+z(n) <= -tau, z free
%      dhigh(i) = - z(i) - w(i)vlow(i), dlow(i) = - z(i) - w(i)vhigh(j), 
%      dhigh(i), dlow(i) free, i = 1,...,n
%      [dlow(i); s(i); a], [dhigh(i);t(i); a] in K_exp, i = 1,...,n

% To write the above problem into the standard form, 
% let the decision variable be 
% [
% [dlow(1); s(i); alow(1)]; ...; [dlow(n); s(n); alow(n)];
% [dhigh(1); s(i); ahigh(1)]; ...; [dhigh(n); s(n); ahigh(n)];
% [w(1); ...; w(n)]; 
% [z(1); ...; z(n)]; zs;
% theta
% ]
% which has total dimension 3*n + 3*n + n + n + 1 + 1 = 8*n + 2

% The additional linear equality constraints are 
% z(1)+...+z(n)+zs = -tau
% ahigh(1) = ahigh(2), ..., ahigh(1) = ahigh(n)
% alow(1) = alow(2), ..., alow(1) = alow(n)
% ahigh(1) = alow(1)

% Generate problem input
n = 50; p = rand(n,1); q = 1-p; vlow = 100 - 50*rand(n,1); vhigh = 100 + 50*rand(n,1);
% Calculate E(V(i))

% save('large_input.mat', 'n', 'p', 'q', 'vlow', 'vhigh');
% load('large_input.mat', 'n', 'p', 'q', 'vlow', 'vhigh'); 
EV = p.*vlow + q.*vhigh; tau = min(110, 0.99*max(EV));

% The original input
% Asset returns
% clear; n = 6; vlow  = [2 -30 -40 -50 -60 -100]'; vhigh  = [2   6   8   10  15   20]'; p = [1 0.055556    0.0625    0.033333    0.053333    0.041667]'; q = 1-p; tau = 3.5;

% Total number of linear equality constraints
m = 5*n + 1;
% RHS vector
b = sparse(m,1); b(end) = -tau; b(end-1) = 1;
% Total dimension of the decision variable
Nt = 3*n + 3*n + n + n + 1;

A = sparse(m, Nt); b = zeros(m,1);
for k = 1:n
    % dlow(i) + z(i) + w(i)*vlow(i) = 0, i = 1,...,n
    row = k; A(row,3*k-2) = 1; A(row,6*n+k) = vlow(k); A(row,7*n+k) = 1;
    % dhigh(i) + z(i) + w(i)*vhigh(i) = 0, i = 1,...,n
    row = n+k; A(row,3*n+3*k-2) = 1; A(row,6*n+k) = vhigh(k); A(row,7*n+k) = 1;
    % p(i)s(i) + q(i)t(i) - alow(1) = 0;
    row = 2*n+k; A(row,3*k-1) = p(k); A(row, 3*n+3*k-1) = q(k); A(row, 3) = -1;
end

for k = 1:n-1
    % alow(1) - alow(i) = 0, i = 2,...,n
    row = 3*n+k; A(row,3) = 1; A(row,3*k+3) = -1;
    % alow(1) - ahigh(i) = 0, i = 2,...,n
    row = 3*n+(n-1)+k; A(row,3) = 1; A(row,3*n+3*k+3) = -1;
end

% alow(1) = ahigh(1)
row = 3*n + 2*(n-1) + 1;
A(row,3) = 1; A(row,3*n+3) = -1;

% w(1)+...w(n) = 1
row = 3*n + 2*(n-1) + 1 + 1;
A(row,3*n+3*n+1:3*n+3*n+n) = ones(1,n); b(row) = 1;

% z(1)+...+z(n)+zs=-tau
row = 3*n + 2*(n-1) + 1 + 1 + 1;
A(row, 3*n+3*n+n+1:3*n+3*n+n+n) = ones(1,n); A(row, 3*n+3*n+n+n+1) = 1; b(row) = -tau;

% Construct cell array input
clear A_cell c_cell blk;
blk{1,1} = 'e'; blk{1,2} = 3*ones(n,1); A_cell{1} = A(:,1:3*n);     c_cell{1} = zeros(3*n,1); c_cell{1}(3) = 1;
blk{2,1} = 'e'; blk{2,2} = 3*ones(n,1); A_cell{2} = A(:,3*n+1:6*n); c_cell{2} = zeros(3*n,1);
blk{3,1} = 'l'; blk{3,2} = n; A_cell{3} = A(:,6*n+1:7*n); c_cell{3} = zeros(n,1);
blk{4,1} = 'u'; blk{4,2} = n; A_cell{4} = A(:,7*n+1:8*n); c_cell{4} = zeros(n,1);
blk{5,1} = 'l'; blk{5,2} = 1; A_cell{5} = A(:,8*n+1);     c_cell{5} = 0;
% blk{6,1} = 'l'; blk{6,2} = 1; A_cell{6} = A(:,8*n+2);     c_cell{6} = 1;

% Call hsd_lqeu
[obj_val, xsol, ysol, zsol, info] = hsd_lqeu(blk, A_cell, c_cell, b, 1e-8, 500);

% Obtain the variables from the stacked optimal solution
a = xsol{1}(3); z = xsol{4}; zs = xsol{5};
dlow = zeros(n,1); dhigh = zeros(n,1); s = zeros(n,1); t = zeros(n,1);
for k = 1:n
    dlow(k) = xsol{1}(3*k-2); s(k) = xsol{1}(3*k-1);
    dhigh(k) = xsol{2}(3*k-2); t(k) = xsol{2}(3*k-1);
end

% Summarize the optimization result
w = xsol{3}; rho_EPSM = 2/(obj_val(1)+obj_val(2));
disp(' '); disp(['The maximum EPSM measure = ' num2str(rho_EPSM) ' with tau = ' num2str(tau) '']);
% disp('The solution (weight) vector = '); disp(w'); 
% disp('The input vhigh, vlow and E(V):'); disp([vhigh'; EV'; vlow']);
disp('Verify that the constraints are satisfied:');
disp([sum(a .* (log(p .* exp(-w .* vlow/a) + q .* exp(-w .* vhigh/a)))) + tau, sum(w)]);
disp('First 10 components of w ='); disp(w(1:min(length(w),10))');