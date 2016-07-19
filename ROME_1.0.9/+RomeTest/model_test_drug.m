% model_test_drug.m
% Drug example from Mosek
% http://www.mosek.com/products/3/tools/doc/html/toolbox/node8.html#SECTION00810000000000000000

% Welcome Message
disp(sprintf('\nTesting Drug Example ... '));

rome_begin;
h = rome_model('Drug Example');

% data 
drug_coeff = [5500, 6100];
raw_coeff  = [-100, -199.9];

% declare variables
drug = rome_model_var(2, 'Cone', rome_constants.NNOC);
raw  = rome_model_var(2, 'Cone', rome_constants.NNOC);

% objective
rome_maximize(raw_coeff*raw + drug_coeff*drug);

% subjected to
p(1) = rome_constraint([0.01*0.995 0.02*0.98]*raw + [-0.5 -0.6]*drug >= 0);
p(2) = rome_constraint(sum(raw) <= 1000);
p(3) = rome_constraint([90 100]*drug <= 2000);
p(4) = rome_constraint([40 50]*drug <= 101);
p(5) = rome_constraint([100 199.9]*raw + [700 800]*drug<=100000);
% p(6) = rome_constraint(drug <= 10);
% p(7) = rome_constraint(raw  <= 200);
% p(8) = rome_constraint(raw  <= 150);

% p(8) = rome_constraint(drug <= 20);

% solve and evaluate
h.solve;
drug_val = h.eval(drug);
raw_val  = h.eval(raw);

% compute difference with old
load +RomeTest/OldDrugData;
disp('Primal Solve... ');
disp(sprintf('Diff in Drug: %g', norm(drug_val - old_drug)));
disp(sprintf('Diff in Raw : %g', norm(raw_val - old_raw)));
disp(sprintf('Diff in Obj : %g', norm(h.objective - old_obj)));

% try solve dual
% h.Solver = 'CPLEXDUAL';
% [x_min, f_min, sol_stat, details] = h.solve;
% disp(sprintf('\nDual Solve... '));
% disp(sprintf('Diff in Drug    : %g', norm(h.eval_var(drug) - old_drug)));
% disp(sprintf('Diff in Raw     : %g', norm(h.eval_var(raw) - old_raw)));
% disp(sprintf('Diff in Obj     : %g', norm(h.objective - old_obj)));
% disp(sprintf('Diff in DualObj : %g', norm(h.DualObjVal - old_objdual)));
% %disp(sprintf('Diff in DualVar : %g', norm(h.eval_var(p) - old_dualvar)));
% disp(sprintf('Diff in DualVar : %g', norm(h.eval_var(p(1:7)) - old_dualvar)));
% % disp(sprintf('DualVars        : %g\n', norm(h.eval_var(p))));


% clear all variables
rome_end;
clearvars -except ROME_ENV


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
