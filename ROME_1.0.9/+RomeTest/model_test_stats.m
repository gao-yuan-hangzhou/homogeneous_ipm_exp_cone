% model_test_stats
% Simple Test Script to test input of statistical properties to uncertain
% variables.

% some parameters
N = 10;
% ind = [2, 4:(N-2)];
ind = 2:N-2;

mean_val = rand(N, 1);
cov_val  = rand(N);
cov_val  = cov_val' * cov_val;
fdev_val = abs(rand(N, 1));
bdev_val = abs(rand(N, 1));
coef_mat = rand(N);

% header
disp(sprintf('\nTesting Variable Statistics ... '));
rome_begin;
h = rome_model('Testing statistics');

% input statistics
z_full = rome_rand_model_var(N);
z_full.set_mean(mean_val);
z_full.Covar = cov_val;
z_full.FDev  = fdev_val;
z_full.BDev  = bdev_val;

% retrieve statistics
dm1 = mean(z_full) - mean_val; dm1 = norm(dm1(:));
ds1 = cov(z_full)  - cov_val;  ds1 = norm(ds1(:));
dp1 = fdev(z_full) - fdev_val; dp1 = norm(dp1(:));
dq1 = bdev(z_full) - bdev_val;  dq1 = norm(dq1(:));
err1 = dm1 + ds1 + dp1 + dq1;
disp(sprintf('Net Error 1 = %0.2f', err1));

% input statistics for a subset
z_subset = rome_rand_model_var(N);
z_subset(ind).set_mean(mean_val(ind));
z_subset(ind).Covar= cov_val(ind, ind);
z_subset(ind).FDev = fdev_val(ind);
z_subset(ind).BDev = bdev_val(ind);

% another way of retrieving statistics
[proj_covar, covar_mix] = h.get_cov(z_subset(ind));
dm2 = mean(z_subset(ind)) - mean_val(ind);      dm2 = norm(dm2(:));
ds2 = proj_covar - covar_mix * cov_val(ind, ind) * covar_mix';  ds2 = norm(ds2(:));
dp2 = fdev(z_subset(ind)) - fdev_val(ind);      dp2 = norm(dp2(:));
dq2 = bdev(z_subset(ind)) - bdev_val(ind);      dq2 = norm(dq2(:));
err2 = dm2 + ds2 + dp2 + dq2;
disp(sprintf('Net Error 2 = %0.2f', err2));

% retrieving statistics on operations
y = coef_mat * z_full;
dm3 = mean(y) - coef_mat * mean_val;            dm3 = norm(dm3(:));
ds3 = 0;
dp3 = 0;
dq3 = 0;
err3= dm3 + ds3 + dp3 + dq3;
disp(sprintf('Net Error 3 = %0.2f', err3));

% clear when done (dont' need to solve)
rome_end;
clearvars -except ROME_ENV;



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
