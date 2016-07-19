function A = solve_sudoku(A_incomplete, b_force_integer)

% +ROMETEST\SOLVE_SUDOKU Helper Routine to solve Sudoku puzzles
%
% Accepts an incomplete Sudoku board as an input and returns a completed 
% Sudoku puzzleboard. Also accepts a scalar flag to decide whether to
% enforce integrality constraints.
%
% Modified by:
% Joel
rome_begin;
h_sudoku_1 = rome_model('Sudoku');
if(nargin == 2 && b_force_integer)
    x = rome_model_var(9, 9, 9, 'Cone', rome_constants.NNOC, 'Continuity', rome_constants.INTEGER);
else
    x = rome_model_var(9, 9, 9, 'Cone', rome_constants.NNOC);
end
rome_constraint(x <= 1);

% Line Constraints
for ii = 1:9
    for jj = 1:9
        rome_constraint(sum(x(ii, jj, :)) == 1);
        rome_constraint(sum(x(ii, :, jj)) == 1);
        rome_constraint(sum(x(:, ii, jj)) == 1);
    end
end

% Box Constraints
for ii = 1:9
    for jj = 1:3
        for kk = 1:3
            y = x(3*(jj-1)+1:3*jj, 3*(kk-1)+1:3*kk, ii);
            rome_constraint(sum(y(:)) == 1);
        end
    end
end

% Data Constraints
for ii = 1:9
    y = x(:, :, ii);
    ind = find(A_incomplete == ii);
    if(isempty(ind)) 
        continue; 
    end
    rome_constraint(y(ind) == 1);
end

rome_minimize(sum(x(:)));
h_sudoku_1.solve;
x_min = h_sudoku_1.eval(x);

A = zeros(9);
for ii = 1:9
    A(find(x_min(:,:,ii))) = ii;
end

rome_end;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
