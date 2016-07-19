% +ROMETEST\MODEL_TEST_SUDOKU Test Script to model and solve Sudoku
% puzzles. 
%
% Calls SOLVE_SUDOKU as core subroutine to solve the actual problem
%
% Modified by:
% Joel

% Create Model Description
disp(sprintf('\nTesting Sudoku Model ... '));

% Create Puzzleboards
A1 = [0 2 0 0 0 0 4 5 0
      3 0 5 0 0 6 0 0 1
      9 0 0 0 8 2 0 0 3
      0 7 0 0 0 0 3 0 5
      0 5 9 0 0 0 8 6 0
      4 0 6 0 0 0 0 9 0
      5 0 0 1 4 0 0 0 7
      8 0 0 6 0 0 1 0 9
      0 9 1 0 0 0 0 4 0];

A2 = [1 0 0 8 0 0 0 3 0
      7 0 0 0 0 0 0 2 0
      0 0 0 5 6 0 0 7 0
      0 0 8 0 0 0 9 0 0
      0 0 5 2 1 7 4 0 0
      0 0 4 0 0 0 7 0 0
      0 3 0 0 8 9 0 0 0
      0 2 0 0 0 0 0 0 8
      0 8 0 0 0 4 0 0 6];

A3 = [0 0 0 0 2 0 9 1 6
      0 0 7 4 8 9 0 0 0
      9 2 5 0 0 0 0 0 8
      0 0 0 1 0 0 4 2 3
      0 9 0 5 3 4 0 7 0
      4 3 1 0 0 2 0 0 0
      8 0 0 0 0 0 2 9 7
      0 0 0 2 6 7 3 0 0
      2 7 3 0 4 0 0 0 0];

A4 = zeros(9);
  
% Call Solver Subroutine to solve
import RomeTest.*;
tic; A1_sol = solve_sudoku(A1, 0); toc; 
tic; A2_sol = solve_sudoku(A2, 0); toc;
tic; A3_sol = solve_sudoku(A3, 0); toc;
tic; A4_sol = solve_sudoku(A4, 1); toc;

% test sudoku against old data
load +RomeTest/OldSudokuData;
disp(sprintf('Error 1 = %g', norm(A1_sol - x_old_sol, 'fro')));
disp(sprintf('Error 2 = %g', norm(A2_sol - x2_old_sol, 'fro')));
disp(sprintf('Error 3 = %g', norm(A3_sol - x3_old_sol, 'fro')));

clearvars -except ROME_ENV;
 
% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
