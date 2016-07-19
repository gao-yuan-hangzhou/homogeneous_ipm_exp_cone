% ROME_FULL_TEST
% SIMPLE TEST SCRIPT

close all;
clear classes; clear global; clear all;

rome_begin;
rome_solver('CPLEX');

% write output to log file
diary(sprintf('./logs/log_%s_%s.txt', datestr(now, 'yyyymmdd'), datestr(now, 'HHMMSS')));
disp(sprintf('\n-------------------------------%%%%%%%%%%%%%%%%%%-------------------------------'));
disp(sprintf('\nExperiment start at %s using %s solver\n', datestr(now), rome_solver));
rome_version;

try
    old_verbose = rome_verbose(0);
    
    % Testing basic operations
    RomeDiagnostics.ops_test;
    
    import RomeTest.*
    model_test_stats;
    model_test_power;
    model_test_beaver_creek;
    model_test_sudoku;
    model_test_norm_min;
    model_test_socp;
    model_test_drug;
    model_test_warehouse;
    model_test_price_of_robustness;
%    model_test_conditional_constraint;
    model_test_seg_inventory;
    model_test_bdldr;
    model_test_project_mgt;
    model_test_inventory;
    
    disp(sprintf('\nExperiment completed at %s\n', datestr(now)));
    disp(sprintf('-------------------------------%%%%%%%%%%%%%%%%%%-------------------------------'));
    diary off;
catch ME
    % turn off diary
    disp(sprintf('\nExperiment Terminated Unsuccessfully at %s\n', datestr(now)));
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!~~~~~~~~~~~~~~~~~~!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    diary off;
    rethrow(ME);
end

rome_verbose(1);

return;
