# Title: Robust Choice
#
# Authors: Shasha Han
# 
#  
#####  Robust Approach on Choice Inconsistency  #####
#
# Test on Cracker dataset, Simulation with |S| = |N|/2
#   
# 
# Shasha Han 
# September, 2015
#
#
#####  SET UP THE PROBLEM  #####
#
param NumObs;
set iniNP := 1..NumObs;  # set of observatios: observation n a is indexed by n
set NP := 1..NumObs/2; # randomly chosen half of the sample
param NALTMAX;
set NALT := 1..NALTMAX;  # set of alternatives
#param NumFixed;
#set NF := 1..NumFixed; # set of explanatory varibles
param gamma; # the maximum number of violations
#param iniGamma;
#set gamma := 1..gammaMAX;
param iniXest{ nn in 1..3292, j in 1..15} >= 0; # the explanatory varialbes
param Xest{ n in NP, j in 1..15} >= 0; # the explanatory varialbes
param ind{ s in 1..1000, nn in iniNP}; # the random permutation of the sample
#param nT;		# number of periods in the data
#set T := 1..nT;	# T is the set of time indices
#param nM;               # number of markets in the data
#set M := 1..nM; # T is the set of markets

var a {n in NP};
var b >=0;
var alpha{j in NALT};
var betaP;
var betaD;
var betaF;

# param iniP {m in M, i in P};

# var p {m in M, i in P} >=0, <=1;

maximize LogLikelihood:  b*( NumObs/2 - gamma) + sum {n in NP} a[n];


subject to PopulationChoices {n in NP, jj in NALT}:
sum { j in NALT } exp ( alpha[j] - alpha[jj] + ( Xest[n, 9+j] - Xest[n, 9+jj] ) * betaP + ( Xest[n, 1+j ] - Xest[n, 1+jj ] ) * betaD + 
                             ( Xest[n, 5+j ] - Xest[n, 5+jj ]) * betaF  + a[n] + (if Xest[n, 14] = jj then 1 else 0 )* b ) <= 1;
Identification: alpha[2] = 0;

####   DEFINE THE PROBLEM   #####

# Name the problem
problem RobustChoice: 

# Choose the objective function
LogLikelihood,

# List the variables
alpha, betaP, betaD, betaF, a, b, 

# List the constraints
PopulationChoices, Identification;


#################################