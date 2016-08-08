addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']); run([fileparts(pwd), '/SDPT3-4.0/startup.m' ]);

clear blk; clear c_cell; clear A_cell; 

blk{1,1} = 'u'; blk{1,2}=2; A_cell{1} = [1,-2]; c_cell{1}=[pi; -2*pi]; b = 1;

% Call hsd_lqeu
display(' '); display('======= Calling hsd_lqeu... ======='); 
hsd_lqeu(blk, A_cell, c_cell, b);
display('======= End calling hsd_lqeu ======='); display(' ');

% Call SDPT3
display(' '); display('======= Calling SDPT3... =======');
sdpt3(blk, A_cell, c_cell, b);
display('======= End calling SDPT3 ======='); display(' ');

% Call sqlp
display(' '); display('======= Calling sqlp... =======');
sqlp(blk, A_cell, c_cell, b);
display('======= End calling sqlp ======='); display(' ');

% Call HSDsqlp
display(' '); display('======= Calling HSDsqlp... =======');
HSDsqlp(blk, A_cell, c_cell, b);
display('End calling HSDsqlp'); display(' ');