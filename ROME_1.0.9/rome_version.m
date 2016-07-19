% ROME_VERSION Displays the version of rome that you are running

function rome_version()

str = sprintf('ROME: Robust Optimization Made Easy\n');
str = [str, sprintf('Version: 1.0.9(beta) \n')];

disp(str);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
