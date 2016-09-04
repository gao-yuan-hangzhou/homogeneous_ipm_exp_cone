clear all

datetime =  datestr(now);
fid = fopen('timestart.txt', 'w');
str = 'the run starts: '
timenow = strcat(str, datetime);
fprintf(fid, '%s', timenow)
fclose(fid); 

% Display welcome message
disp('Robust Choice');
% Modified by Shasha, Aug 22st, 2015

% Declare GLOBAL variables
% GLOBAL variables are all in caps
% DO NOT CHANGE ANY OF THESE 'global' STATEMENTS
global NP NALTMAX  ScenarioCase 
global PMAT DMAT FMAT CMAT

delete('Sample_estimates.txt');

disp('Importing data arrays for run.');

load SamplePrice; % load the price data 3291 * 4
load SampleDisp; % load the disp 3292 * 4
load SampleFeat;
load SampleChoice7;
PMAT = SamplePrice; 
DMAT = SampleDisp;
FMAT = SampleFeat;
CMAT = SampleChoice7;


ScenarioCase = 5;
NP = 100; % Number of people (decision-makers) in dataset 
NALTMAX = 4; % Number of alternatives in the dataset


d = [0.15, 0.15, 0.45, 0.25];
for sn = 0 : ScenarioCase 
    N = sum( CMAT( (sn * NP + 1) : (( sn + 1) * NP), : ), 1);
    gamma = floor (N .* (1 - d));  
% Set up parameter estimation problem
h = rome_begin('Parameter Estimation Under Robust Choice');

newvar alpha(1,NALTMAX)
newvar s(NP,NALTMAX)
newvar t(NP,NALTMAX)
newvar a(NP)
newvar beta(3)
newvar b(1,NALTMAX)


 rome_minimize( -sum(a) - sum (b .* ( N - gamma)))
 rome_constraint (b >= 0);
 
for  n = (sn * NP + 1) : (( sn + 1) * NP)
     m = n - (sn * NP);
     
     alpha(2) = 0;
     alpha(1) = 0 ;
      alpha(3) = 0;
     alpha(4) = 0 ;
     
     beta(2) = 0;
     beta(3) = 0;
      
      for j =1:NALTMAX
          if j == j* CMAT(n,j)
              for i =1:NALTMAX
              rome_constraint(t(m,i) >=  exp( alpha(i)- alpha(j) + ( PMAT(n,i)- PMAT(n, j) ) * beta(1) + ( DMAT(n,i)- DMAT(n, j) ) * beta(2) + ( FMAT(n,i)- FMAT(n, j) ) * beta(3) + a(m) + b (j) )); % here beta is the coffefficent
              end
              rome_constraint(sum(t(m,:)) <= 1);
          else 
              for i =1:NALTMAX
              rome_constraint(s(m,i) >=  exp( alpha(i)- alpha(j) + ( PMAT(n,i)- PMAT(n, j) ) * beta(1) + ( DMAT(n,i)- DMAT(n, j) ) * beta(2) + ( FMAT(n,i)- FMAT(n, j) ) * beta(3) + a(m) ));
              end
              rome_constraint(sum(s(m,:)) <= 1);
          end
         
      end
    
end

 
      
h.solve('CPLEX'); 
beta = h.eval(beta);  % Output solution
alpha = h.eval(alpha);

NegLogLihood  = h.objective; % Output negative loglikelihood

rome_end;
%disp('Estimatation Parameter');
%beta
%alpha
%disp('Negative loglikelihood');
%NegLogLihood

%%%%% write the estimates in each subsample
fidx = fopen('Sample_estimates.txt', 'at');
fprintf(fidx, '%g %g %g %g %g %g %g %g \n', alpha, beta, NegLogLihood);
fprintf('\n');
fclose(fidx);

end   
datetime =  datestr(now);
fid = fopen('timefinish.txt', 'w');
str = 'the run finishs :'
timenow = strcat(str, datetime);
fprintf(fid, '%s', timenow)
fclose(fid); 