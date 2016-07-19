% ops_test
% Script to test correctness of basic operations: i.e. subsituting some
% values and making sure that at the end of the day everything corresponds
% nicely. 

h1 = rome_model('Test Basic Operations');

% output data structure
test_array = [];
descr_array = {};

% perform tests
import RomeDiagnostics.*;
tic; ops_test_mtimes; toc;      % 1
tic; ops_test_times; toc;       % 2
tic; ops_test_rdivide; toc;     % 3
tic; ops_test_plus; toc;        % 4
tic; ops_test_minus; toc;       % 5
tic; ops_test_transpose; toc;   % 6
tic; ops_test_sum; toc;         % 7
tic; ops_test_rand_var; toc;    % 8
tic; ops_test_cat; toc;         % 9
tic; ops_test_assign; toc;      % 10
tic; ops_test_delete; toc;      % 11

disp(sprintf('Tests Complete!\n'));

% report test results
if(length(test_array) ~= length(descr_array))
    error('Number of tests Mismatch!');
end
report_str = '';
for ii = 1:length(test_array)
    report_str = [report_str, sprintf('%s    = %0.12E\n', descr_array{ii}, test_array(ii))];
end
disp(report_str);
disp(sprintf('Net Error: %f', norm(test_array)));

clearvars -except h1 ROME_ENV

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
