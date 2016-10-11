addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;

global N Y J Q disp_mat feat_mat price_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set (total) sample size, Gamma/N and partition ratio
sample_size = 200;
% This means use(1-partition_ratio) data for training and partition_ratio for testing
partition_ratio = 0.5;
% This means we vary Gamma according to the list
% Useless since we are solving the case with Gamma=0
Gamma_list = (0:10);
Gamma = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basic model setting: there is a customer and J types of products. 
% Each time the customer purchases one of the products.
% There have been N purchases.

% Assign codes to the brands
code_sunshine = 1; code_keebler = 2;
code_nabisco = 3; code_private = 4;

% J is the number of product types
J = 4;

% Read the prucahse record data into a matrix
purchase_record_mat = csvread('data/Cracker_matlab_friendly.csv',2);
total_num_obs = size(purchase_record_mat,1);
purchase_record_mat = purchase_record_mat(randsample(total_num_obs,min(total_num_obs,sample_size)),:); % Get first 10 observations
%save('purchase_record_mat.mat', 'purchase_record_mat');
load('purchase_record_mat.mat');

% N is the number of observations of this customer
Ntotal = size(purchase_record_mat,1);
disp(['Total number of observations = ' num2str(Ntotal)]);

% Q is the number of attributes
Q = 3;

idx = 1;
% Partition the data into training and testing data
purchase_record_mat_train = purchase_record_mat(1:floor(sample_size*(1-partition_ratio)),:);
purchase_record_mat_test = purchase_record_mat(floor(sample_size*(1-partition_ratio))+1:end,:);

% Construct the N-by-J matrix representing this customer's purchase reord
Y = purchase_record_mat_train(:,end-J+1:end);
N = size(purchase_record_mat_train,1);
disp(['Number of observations in the training set = ' num2str(N)]);

% Compute N-dimensional j_hat, which is defines as follows
% j_hat(n) = k means in the n-th observation the customer bought item k, n=1:N, k=1:J.
j_hat = Y*(1:J)';
% Construct N-by-J matrices representing "disp", "feat" and "price" attributes
% In the model, the vector x(n,j) = [disp_mat(n,j); feat_mat(n,j); prie_mat(n,j)]
disp_mat = purchase_record_mat_train(:,2:1+J);
feat_mat = purchase_record_mat_train(:,2+J:1+2*J);
price_mat = purchase_record_mat_train(:,2+2*J:1+3*J);


f_min_vec = @(abvec) f_min(abvec(1:J), abvec(J+1:J+Q));

[xopt,fval] = fminunc(f_min_vec,zeros(J+Q,1));
alpha_opt1 = xopt(1:J);
beta_opt1 = xopt(J+1:J+Q);




