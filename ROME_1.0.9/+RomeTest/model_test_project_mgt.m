% +ROMETEST\MODEL_TEST_PROJECT_MGT Test Script to model and solve
% the project management problem.
%
% Calls the following subroutines to solve the actual problem:
%   1. SOLVE_PROJECT_MGT_LDR (solve using linear decision rule)
%   2. SOLVE_PROJECT_MGT_DLDR (solve using deflected linear decision rule)
%   3. SOLVE_PROJECT_MGT_SLDR (solve using segregated linear decision rule)
%   4. SOLVE_PROJECT_MGT_SDLDR (solve using segregated-deflected linear decision rule)
%
% References: 
% 1.  X. Chen, M. Sim, P. Sun, and J. Zhang, "A Linear-Decision Based 
%     Approximation Approach to Stochastic Programming", Operations 
%     Research, 2008, 56(2), 344-357.
%
% Modified by:
% 1. Joel (30 Oct 2008)
%

% Import calling routines
import RomeTest.*;

% display welcome
disp(sprintf('\nTesting Project Management Example...'));

C_arr = [8, 19];      % project budget
beta_arr = [0.1, 0.2, 0.3, 0.4]; % parameter used to define uncertainty support

W = 6;  % Width of project grid
H = 4;  % Height of project grid
a = 3;  % linear part for relating time with resources
b = 3;  % constant in relating time with resources
c = 1;  % cost of each unit of resource (constant)

% Parameter for Segregated Decision Rules
z_pos_mean = 0.5;   % mean of positive part of primitive uncertainty

% result arrays
result_mat_ldr = zeros(numel(C_arr), numel(beta_arr));
result_mat_dldr = zeros(numel(C_arr), numel(beta_arr));
result_mat_dldr_auto = zeros(numel(C_arr), numel(beta_arr));
result_mat_sldr = zeros(numel(C_arr), numel(beta_arr));
result_mat_sdldr = zeros(numel(C_arr), numel(beta_arr));
result_mat_sdldr_auto = zeros(numel(C_arr), numel(beta_arr));

runtime_ldr = zeros(numel(C_arr), numel(beta_arr));
runtime_dldr = zeros(numel(C_arr), numel(beta_arr));
runtime_dldr_auto = zeros(numel(C_arr), numel(beta_arr));
runtime_sldr = zeros(numel(C_arr), numel(beta_arr));
runtime_sdldr = zeros(numel(C_arr), numel(beta_arr));
runtime_sdldr_auto = zeros(numel(C_arr), numel(beta_arr));

% override MOSEK PARAM for stability
rome_begin;
global ROME_ENV;
ROME_ENV.MSK_PARAMS.MSK_DPAR_INTPNT_TOL_PSAFE = 1E9;

% Iterate over each C and each beta
for ii = 1:numel(C_arr)
    for jj = 1:numel(beta_arr)
        
        % display progress
        disp(sprintf('Iter (%d, %d): C = %0.2f, beta = %0.2f', ii, jj, C_arr(ii), beta_arr(jj)));
        
        % RUN Linear Decision rule
        tic;
        [min_time_ldr, x_val_ldr, y_val_ldr] = solve_project_mgt_ldr(C_arr(ii), beta_arr(jj), W, H, a, b, c);
        runtime_ldr(ii, jj) = toc;
        result_mat_ldr(ii, jj) = min_time_ldr;

        % RUN Deflected Linear Decision rule
        tic;
        [min_time_dldr, x_val_dldr, y_val_dldr, r_val_dldr] = solve_project_mgt_dldr(C_arr(ii), beta_arr(jj), W, H, a, b, c);
        runtime_dldr(ii, jj) = toc;
        result_mat_dldr(ii, jj) = min_time_dldr;
        
        % RUN Auto Deflected Linear Decision rule
        tic;
        [min_time_dldr_auto, x_val_dldr_auto, y_val_dldr_auto] = ...
            solve_project_mgt_dldr_auto(C_arr(ii), beta_arr(jj), W, H, a, b, c);
        runtime_dldr_auto(ii, jj) = toc;
        result_mat_dldr_auto(ii, jj) = min_time_dldr_auto;

        if(strcmp(rome_solver(), 'SDPT3DUAL'))
           continue; 
        end
        
        % RUN Segregated Linear Decision rule
        tic;
        [min_time_sldr, x_val_sldr, y_val_sldr] = solve_project_mgt_sldr(C_arr(ii), beta_arr(jj), W, H, a, b, c, z_pos_mean);
        runtime_sldr(ii, jj) = toc;
        result_mat_sldr(ii, jj) = min_time_sldr;
 
        % RUN Segregated-Deflected Linear Decision rule
        tic;
        [min_time_sdldr, x_val_sdldr, y_val_sdldr] = solve_project_mgt_sdldr(C_arr(ii), beta_arr(jj), W, H, a, b, c, z_pos_mean);
        runtime_sdldr(ii, jj) = toc;
        result_mat_sdldr(ii, jj) = min_time_sdldr;
 
        % RUN AUTO Segregated-Deflected Linear Decision rule
        tic;
        [min_time_sdldr_auto, x_val_sdldr_auto, y_val_sdldr_auto] = ...
                solve_project_mgt_sdldr_auto(C_arr(ii), beta_arr(jj), W, H, a, b, c, z_pos_mean);
        runtime_sdldr_auto(ii, jj) = toc;
        result_mat_sdldr_auto(ii, jj) = min_time_sdldr_auto;
    end
end

% override MOSEK PARAM for stability
% restore default
ROME_ENV.MSK_PARAMS = [];

% check vs benchmark
load +RomeTest/ProjMgtData;

% number of iterations
Niter = numel(C_arr) * numel(beta_arr);

% display errors
cur_err = result_mat_ldr - result_ldr_benchmark;
disp(sprintf('     LDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err(:)), sum(runtime_ldr(:)) ./ Niter));

cur_err_dldr = result_mat_dldr - result_dldr_benchmark;
disp(sprintf('    DLDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_dldr(:)), sum(runtime_dldr(:)) ./ Niter));

cur_err_dldr_auto = result_mat_dldr_auto - result_dldr_benchmark;
disp(sprintf(' DLDR(A): Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_dldr_auto(:)), sum(runtime_dldr_auto(:)) ./ Niter));

cur_err_sldr = result_mat_sldr - result_sldr_benchmark;
disp(sprintf('    SLDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_sldr(:)), sum(runtime_sldr(:)) ./ Niter));

cur_err_sdldr = result_mat_sdldr - result_sdldr_benchmark;
disp(sprintf('   SDLDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_sdldr(:)), sum(runtime_sdldr(:)) ./ Niter ));

cur_err_sdldr_auto = result_mat_sdldr_auto - result_sdldr_benchmark;
disp(sprintf('SDLDR(A): Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_sdldr_auto(:)), sum(runtime_sdldr_auto(:)) ./ Niter ));

% Build comparison table
C_tmp      = repmat(C_arr, numel(beta_arr), 1);
beta_tmp   = repmat(beta_arr, 1, numel(C_arr));
result_tmp = result_mat_ldr';
result_dldr_tmp = result_mat_dldr';
result_dldr_auto_tmp = result_mat_dldr_auto';
result_sldr_tmp = result_mat_sldr';
result_sdldr_tmp = result_mat_sdldr';
result_sdldr_auto_tmp = result_mat_sdldr_auto';

% Make table
T = [C_tmp(:), beta_tmp(:), result_tmp(:), result_dldr_tmp(:), result_dldr_auto_tmp(:), ...
    result_sldr_tmp(:), result_sdldr_tmp(:), result_sdldr_auto_tmp(:)];

% Make output string
str = sprintf('%8s   |%8s   |%8s   |%8s   | %8s  |%8s   |%8s   |  %8s\n', ...
               'C', 'beta', 'Z_LDR', 'Z_DLDR', 'Z_DLDRA', 'Z_SLDR', 'Z_SDLDR', 'Z_SDLDRA'); 
str = [str, repmat('-', 1, length(str)), sprintf('\n')]; 
for ii = 1:size(T, 1)
    str = [str, sprintf('%8.2f   |%8.3f   |%8.3f   |%8.3f   |%8.3f   |%8.3f   |%8.3f   |%8.3f\n',...
                        T(ii, 1), T(ii, 2), T(ii, 3), T(ii, 4), T(ii, 5), T(ii, 6), T(ii, 7), T(ii, 8))]; 
end
disp(str);

% clear all data when done
clearvars -except ROME_ENV


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
