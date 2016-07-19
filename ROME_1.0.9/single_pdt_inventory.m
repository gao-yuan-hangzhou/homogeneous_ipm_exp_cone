% ROME Example: Single Product Multi-Period Inventory Control
%
% See: Single Product Multi-Period Inventory Control.pdf
%
% Modified by:
% 1. Joel & Melvyn   (Created 27 July 2009)

% Display welcome message
disp('Single Product Multi-Period Inventory Control');

% Inventory setup parameters
T = 10;                     % Number of periods 
xMax = 260;                 % Maximum ordering level
c = 0.1*ones(T,1);          % Ordering costs
h = 0.02*ones(T,1);         % Holding costs
b = [0.2*ones(T-1,1); 2];   % Shortage costs

% Demand information
mu = 200;       % mean 
alpha = 0.5;    % degree of demand correlation
zRange = mu/T;
zSample = -zRange:zRange/20:zRange;           
zDist   = ones(length(zSample),1)/length(zSample);
% Determine deviation information using divmeasure function
[zMean, zStd, zFDev, zBDev,zMin,zMax] =  divmeasure(zSample,zDist);

% begin rome
rome_h = rome_begin;   

% Define uncertain factors
newvar z(T) uncertain;
z.set_mean(0);
rome_box(z, zMin, zMax);
z.Covar = zStd^2;  % Covariance is diagonal with values equal to variance
z.FDev = zFDev;
z.BDev = zBDev;

% Demand relations with uncertain factors
d = mu + alpha*z(1);
for t = 2:T
    d(t) = d(t-1) - (1-alpha)*z(t-1) + z(t);
end

% Model variables
% Ordering levels
pX = logical([tril(ones(T)), zeros(T, 1)]); % Dependency patter of x on z. 
newvar x(T, z, 'Pattern', pX) linearrule;

% Inventory level
newvar y(T+1,z) linearrule;   % Note that dependency of y is enforced in
                              % inventory balance equation
% Objective
rome_minimize(c'*mean(x) + h'*mean(pos(y(2:T+1))) + b'*mean(neg(y(2:T+1))));

% Constraints
rome_constraint(y(1)==0);  % Initial inventory level
for t = 1:T
    rome_constraint(y(t+1)==y(t) + x(t) - d(t));
end

% order quantity constraint
rome_box(x, 0, xMax);

% solve
rome_h.solve_deflected;
obj_val = rome_h.objective;
x_val   = rome_h.eval(x);

rome_end;

disp(obj_val);
x_val


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
