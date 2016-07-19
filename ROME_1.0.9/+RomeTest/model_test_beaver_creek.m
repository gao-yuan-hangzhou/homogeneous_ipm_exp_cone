% model_test_beaver_creek
% Simple Test Program for 2 variable LP.

% Welcome message
disp(sprintf('\nTesting Beaver Creek Model ...'));

% Expt 1: Simple Beaver Creek Problem from Taylor Textbook
rome_begin;
tic;
h_beaver_creek = rome_model('Beaver Creek');

% set up modeling variables
x = rome_model_var('Cone', rome_constants.NNOC);
y = rome_model_var('Cone', rome_constants.NNOC);

% set up objective function
rome_maximize(40 * x +   50 * y);

% constraints
p = rome_constraint(x + 2*y <= 40);
q = rome_constraint(4*x + 3*y <= 120);

% solve
h_beaver_creek.solve;
t = toc;

% report solution
disp(sprintf('Solution = (%g, %g), Obj = %g, time taken = %g', ...
        h_beaver_creek.eval(x), h_beaver_creek.eval(y), h_beaver_creek.objective, t));
    
% Expt 2: changing constraint RHS
tic;
pRHS = linspace(35, 10, 20);
disp(sprintf('\nBegin Changing Constraint 1 RHS from %0.2f to %0.2f: ', pRHS(1), pRHS(end)));
for ii = 1:length(pRHS)
    % modify the constraint
    rome_mod_constraint(p, x + 2*y <= pRHS(ii))
    
    % re-solve
    h_beaver_creek.solve;

    % report solution
    disp(sprintf('Iter %0.2d, pRHS = %0.4f: Solution = (%0.3f, %0.3f), Obj = %0.2f', ...
        ii, pRHS(ii), h_beaver_creek.eval(x), h_beaver_creek.eval(y), h_beaver_creek.objective));
end
% undo previous change
rome_mod_constraint(p, x + 2*y <= 40); 
t = toc;
disp(sprintf('All Complete, time taken = %g', t));
    
% Expt 3: changing constraint coeff
tic;
qCoeff = linspace(4, 10, 25);
disp(sprintf('\nBegin Changing Constraint 2 x-coefficient from %0.2f to %0.2f: ', qCoeff(1), qCoeff(end)));
for ii = 1:length(qCoeff)
    % modify constraint
    rome_mod_constraint(q, qCoeff(ii) * x + 3*y <= 120);
    
    % re-solve
    h_beaver_creek.solve;

    % report solution
    disp(sprintf('Iter %0.2d, qXCoeff = %0.4f: Solution = (%0.3f, %0.3f), Obj = %0.2f', ...
        ii, qCoeff(ii), h_beaver_creek.eval(x), h_beaver_creek.eval(y), h_beaver_creek.objective));
end
t = toc;
disp(sprintf('All Complete, time taken = %g', t));

% Expt 4: changing both constraint coeff
% undo previous change
rome_mod_constraint(q, 4*x + 3*y <= 120); 
tic;
qCoeff = linspace(4, 10, 5);
pCoeff = linspace(3, 15, 5);

min_xy  = NaN;
min_qp  = NaN;
min_romeit = Inf;

disp(sprintf('\nBegin Changing Constraint 2 x-coefficient from %0.2f to %0.2f and Constraint 2 y-coefficient from %0.2f to %0.2f: ', ...
    qCoeff(1), qCoeff(end), pCoeff(1), pCoeff(end)));

for ii = 1:length(qCoeff)
    for jj = 1:length(pCoeff)
        % modify constraint
        rome_mod_constraint(q, qCoeff(ii) * x + pCoeff(jj)*y <= 120);

        % re-solve
        h_beaver_creek.solve;

        % report solution
        disp(sprintf('Iter %0.2d, qXCoeff = %0.4f, pYCoeff = %0.4f: Solution = (%0.3f, %0.3f), Obj = %0.2f', ...
            (ii-1) * length(pCoeff) + jj , qCoeff(ii), pCoeff(jj), h_beaver_creek.eval(x), h_beaver_creek.eval(y), h_beaver_creek.objective));
        
        % store worst-case results
        if(h_beaver_creek.objective < min_romeit)
            min_romeit = h_beaver_creek.objective;
            min_xy = [h_beaver_creek.eval(x), h_beaver_creek.eval(y)];
            min_qp = [qCoeff(ii), pCoeff(jj)];
        end
    end
end

disp(sprintf('\nSummary: \nWorst-case(qXCoeff, pYCoeff) = (%0.4f, %0.4f) \nWorst-case Sol (X, Y) = (%0.3f, %0.3f) \nWorst-case Objective = %0.2f', ...
     min_qp(1), min_qp(2), min_xy(1), min_xy(2), min_romeit));
t = toc;
disp(sprintf('All Complete, time taken = %g', t));


% Expt 5: using robust optimization
% undo previous change
rome_mod_constraint(q, 4*x + 3*y <= 120); 
tic;

% add coefficients and constraints
cx = rome_rand_model_var;
rome_constraint(cx >= 4);
rome_constraint(cx <= 10);

cy = rome_rand_model_var;
rome_constraint(cy >= 3);
rome_constraint(cy <= 15);

% add robust constraint
rome_constraint(cx * x + cy * y <= 120);

% solve
h_beaver_creek.solve;
t = toc;
% report solution
disp(sprintf('\nRobust Solution = (%g, %g), Obj = %g, time taken = %g', ...
        h_beaver_creek.eval(x), h_beaver_creek.eval(y), h_beaver_creek.objective, t));

rome_end;
% clear vars
clearvars -except ROME_ENV;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
