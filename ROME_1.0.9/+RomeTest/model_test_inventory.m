% +ROME_EXAMPLES\MODEL_TEST_INVENTORY
%
% Script to test various values of the inventory problem under various
% decision rules. To ensure that values agree with those in the paper by
% See and Sim (2009).
%
% Modified by:
% 1. Joel (Created 20 May 2009)

% display welcome message
disp('Testing Inventory Problem ... ');

c0 = 0.1;
h0 = 0.02;
xMax = 260;

alpha_vec = 0:0.25:1;
boverh_vec = 10:20:50;
N = 5;

result_tldr = zeros(length(alpha_vec), length(boverh_vec));
time_tldr   = zeros(length(alpha_vec), length(boverh_vec));

result_dldr = zeros(length(alpha_vec), length(boverh_vec));
time_dldr   = zeros(length(alpha_vec), length(boverh_vec));

result_ldr = zeros(length(alpha_vec), length(boverh_vec));
time_ldr   = zeros(length(alpha_vec), length(boverh_vec));

result_sdr = zeros(length(alpha_vec), length(boverh_vec));
time_sdr   = zeros(length(alpha_vec), length(boverh_vec));

import RomeTest.*;

for ii = 1:size(result_tldr, 1)
    for jj = 1:size(result_tldr, 2)
        % instantiate parameters
        if(N <= 10)
            mu = 200;
        else
            mu = 240;
        end
        hcost = h0 * ones(N, 1);
        pcost = boverh_vec(jj) * hcost;
        pcost(end) = 10 * pcost(1);
        alpha = alpha_vec(ii);
        c = c0 * ones(N, 1);
        
        % display
        fprintf('Iter (%d, %d): alpha = %0.2f, b/h = %2g \n', ii, jj, alpha, boverh_vec(jj));
        
        % solve Truncated Linear Decision Rule
        tic;
        result_tldr(ii, jj) = solve_robust_inventory_tldr(N, c, hcost, pcost, mu, xMax, alpha);
        time_tldr(ii, jj) = toc;
        
        % solve Deflected Linear Decision Rule
        tic;
        result_dldr(ii, jj) = solve_robust_inventory_dldr(N, c, hcost, pcost, mu, xMax, alpha);
        time_dldr(ii, jj) = toc;
        
        % solve Linear Decision Rule
        tic;
        result_ldr(ii, jj) = solve_robust_inventory_ldr(N, c, hcost, pcost, mu, xMax, alpha);
        time_ldr(ii, jj) = toc;
        
        % solve Static Decision Rule
        tic;
        result_sdr(ii, jj) = solve_robust_inventory_sdr(N, c, hcost, pcost, mu, xMax, alpha);
        time_sdr(ii, jj) = toc;
    end
end

% compare against benchmark
load +RomeTest/InventoryData.mat;

% display errors
Niter = numel(result_tldr);
cur_err_tldr = result_tldr - result_tldr_benchmark;
disp(sprintf('TLDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_tldr(:)), sum(time_tldr(:)) ./ Niter));

cur_err_dldr = result_dldr - result_dldr_benchmark;
disp(sprintf('DLDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_dldr(:)), sum(time_dldr(:)) ./ Niter));

cur_err_ldr = result_ldr - result_ldr_benchmark;
disp(sprintf(' LDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_ldr(:)), sum(time_ldr(:)) ./ Niter));

cur_err_sdr = result_sdr - result_sdr_benchmark;
disp(sprintf(' SDR: Error = %0.3f, Avg Time Taken: %0.3f secs.', norm(cur_err_sdr(:)), sum(time_sdr(:)) ./ Niter));

% Make table
alpha_tmp  = repmat(alpha_vec, numel(boverh_vec), 1);
boverh_tmp = repmat(boverh_vec, 1, numel(alpha_vec));
result_tldr = result_tldr';
result_dldr = result_dldr';
result_ldr = result_ldr';
result_sdr = result_sdr';


T = [alpha_tmp(:), boverh_tmp(:), result_tldr(:), result_dldr(:), result_ldr(:), result_sdr(:)];

% Make output string
str = sprintf('%8s   |%8s   |%8s   |%8s   |%8s   |%8s   \n', ...
               'alpha', 'b/h', 'Z_TLDR', 'Z_DLDR', 'Z_LDR', 'Z_SDR'); 
str = [str, repmat('-', 1, length(str)), sprintf('\n')]; 
for ii = 1:size(T, 1)
    str = [str, sprintf('%8.2f   |%8.3f   |%8.3f   |%8.3f   |%8.3f   |%8.3f\n', ...
                        T(ii, 1), T(ii, 2), T(ii, 3), T(ii, 4), T(ii, 5), T(ii, 6))]; 
end
disp(str);

% clear all data when done
clearvars -except ROME_ENV


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
