% ROMETEST\MODEL_TEST_SIMPLE_LP Script file solving a simple Linear
% Program.

% Set up ROME Environment
rome_begin; 
h = rome_model('Simple LP'); % Create Rome Model

% set up modeling variables
newvar x y;

% set up objective function
rome_maximize(12 * x +  15 * y);

% input constraints constraints
rome_constraint(x + 2*y <= 40);
rome_constraint(4*x + 3*y <= 120);
rome_constraint(x >= 0);
rome_constraint(y >= 0);

% solve
h.solve;

% extract objective values
x_val = h.eval(x);
y_val = h.eval(y);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
