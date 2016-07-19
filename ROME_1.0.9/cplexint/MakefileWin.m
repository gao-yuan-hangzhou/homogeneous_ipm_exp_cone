%=========================
% CPLEX VERSION 10.0 !!!!
%=========================
disp('Automatic compilation of CPLEXINT: a Matlab MEX interface for CPLEX90 under Windows.');
disp(' ');
disp('We assume that CPLEX90 is installed in: C:\ILOG\cplex90');
disp('and that you are using MSVC6 as a compiler in Matlab (see mex -setup).');
disp('If that is not correct path for CPLEX90, modify this file accordingly.');
disp(' ');
disp('***************************************************************');
disp('*  Old version of MEX CPLEXINT files will be overwritten !!!  *');
disp('***************************************************************');
disp(' ');
disp('Press any key to continue');
disp(' ');
pause;
disp('Compiling in progress. Please wait.');

% Windows environment
%--------------------

% compile cplexin that uses cplex90
mex -I'C:\ILOG\cplex112\include\ilcplex' cplexint.c 'C:\ILOG\CPLEX112\lib\x86_.net2005_8.0\stat_mda\cplex112.lib'
%mex -win64 -I'C:\ILOG\CPLEX122\cplex\include\ilcplex' cplexint.c 'C:\ILOG\CPLEX122\cplex\lib\x86_.net2005_8.0\stat_mda\cplex122.lib'

disp('Compiling finished.');

return;
