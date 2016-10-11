%%*************************************************************************
%% Run this script in Matlab command window to generate some mex files
%% 
%%*************************************************************************
%% Copyright (c) 2014 by
%% Kim-Chuan Toh
%%*************************************************************************

   function Installmex 

   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4));
%%
   if strcmp(computer_model,'PCWIN')
      str0 = ['''',matlabroot,'\extern\lib\win32\lcc\'' '];
      if (exist(eval(str0),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
      else
         str1 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwblas.lib''  '];
      end
      libstr = [str1,str2];
   elseif strcmp(computer_model,'PCWIN64')
      str0 = ['''',matlabroot,'\extern\lib\win64\lcc\'' '];
      if (exist(eval(str0),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwblas.lib''  '];
      else
         str1 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
      end
      libstr = [str1,str2];
   else
      libstr = '  -lmwlapack -lmwblas  '; 
   end
   mexcmd = 'mex -O  -largeArrayDims  -output ';    
%%
   if (matlabversion < 7.3)  && (matlabrelease <= 2008)
      error(' needs MATLAB version 7.4 and above'); 
   end
   fsp = filesep;

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir); 
%%
%% generate mex files in mexfun
%%
   clear fname

   src = [curdir,fsp,'Mexfiles']; 
   eval(['cd ','Mexfiles']); 
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);       
   %%
   fname{1} = 'mexbwsolve';
   fname{2} = 'mexfwsolve';
   fname{3} = 'mextriang'; 
   fname{4} = 'mexMatvec';
   fname{5} = 'mexFnorm';
 
   for k = 1:length(fname)
       cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]);  
   end         
    
   fprintf ('\n Compilation of mexFunctions completed.\n'); 
   cd .. 
  
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
