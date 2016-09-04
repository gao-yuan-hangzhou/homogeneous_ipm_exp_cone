addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;

% Basic model setting: there is a customer and J types of products. 
% Each time the customer purchases one of the products.
% There have been N purchases.

% Assign codes to the brands
code_sunshine = 1; code_keebler = 2;
code_nabisco = 3; code_private = 4;

% J is the number of product types
J = 4;

% Consider ONE customer with a specific id
curr_id = 33;

% Read the prucahse record data into a matrix
% Each row of the matrix is one record in the following format:
% 
purchase_record_mat = csvread('data/Cracker_matlab_friendly.csv',2);

% Get the submatrix of this customer's purcahse record
id_repeat = purchase_record_mat(:,1);
indices_curr_id = find(id_repeat==curr_id);
curr_id_record_mat = purchase_record_mat(indices_curr_id,:);

% N is the number of observations of this customer
N = length(indices_curr_id);

% Now construct the N-by-J matrix representing this customer's purchase reord
Y = curr_id_record_mat(:,end-J+1:end);

% Compute N-dimensional j_hat, which is defines as follows
% j_hat(n) = k means in the n-th observation the customer bought item k, n=1:N, k=1:J.
j_hat = Y*(1:J)';

% Construct N-by-J matrices representing "disp", "feat" and "price" attributes
% In the model, the vector x(n,j) = [disp_mat(n,j); feat_mat(n,j); prie_mat(n,j)]
disp_mat = curr_id_record_mat(:,2:1+J);
feat_mat = curr_id_record_mat(:,2+J:1+2*J);
price_mat = curr_id_record_mat(:,2+2*J:1+3*J);

