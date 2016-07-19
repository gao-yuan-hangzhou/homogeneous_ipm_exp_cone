% ROME Example: Portfolio Allocation under the Entropic Satisficing Measure
%
% Modified by:
% 1. Melvyn   (Created 17 Aug 2009)

clear
% Display welcome message
disp('Portfolio Allocation under the Entropic Satisficing Measure');

% Asset returns
N = 6;  % Number of assets
% Asset returns in percentage
vlo  = [2  -30 -40 -50 -60 -100];
vhi  = [2   6   8   10  15   20];
% probability of asset lowest value
p = [1 0.055556    0.0625    0.033333    0.053333    0.041667];
q = 1-p;

tau = 3.5; %Target return


% Set up portfolio optimization problem
h = rome_begin('Portfolio Allocation under the Entropic Satisficing Measure');
newvar('w(N)', 'nonneg');
newvar('a', 'nonneg');
newvar z(N) 

rome_minimize(a)
rome_constraint(sum(z)<=-tau);

for i=1:N
   rome_constraint(a >= p(i)*expcone( -z(i)-vlo(i)*w(i),a) + q(i)*expcone(-z(i)-vhi(i)*w(i),a));
end;
rome_constraint(sum(w)==1);
 
h.solve('SDPT3'); 
wsol = h.eval(w);  % Output solution
rho  = 1/h.objective; % Output EPSM value

rome_end;

disp(['Target=' num2str(tau)]);
disp(['EPSM level = ' num2str(rho)]);
disp('Solution (weight) vector:');
disp(wsol');

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.