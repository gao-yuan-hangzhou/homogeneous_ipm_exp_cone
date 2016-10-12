addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set (total) sample size, Gamma/N and partition ratio
sample_size = 100;
% This means use(1-partition_ratio) data for training and partition_ratio for testing
partition_ratio = 0.5;
% This means we vary Gamma according to the list
Gamma_list = (0:5); 
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

% For hsd_lqeu input, the decision variables are
% alpha, beta, a[n], b>=0, [s_njl; t_njl; u_njl] in K_exp
dim_x = (J+Q) + N + 1 + 3*N*J^2;

% Construct the cell array blk
blk{1,1} = 'u'; blk{1,2} = J; % alpha
blk{2,1} = 'u'; blk{2,2} = Q; % beta
blk{3,1} = 'u'; blk{3,2} = N; % dual variable a[n], n = 1:N
blk{4,1} = 'l'; blk{4,2} = 1; % dual variable b
blk{5,1} = 'e'; blk{5,2} = 3*ones(N*J^2,1); % aux exp cone variables [s,t,u]

% Construct the cell array c
c{1} = zeros(J,1);
c{2} = zeros(Q,1);
c{3} = -ones(N,1);
%c{4} = -(N - Gamma);
c{5} = zeros(3*(N*J^2),1);
% Now assemble part of the matrix A in Ax=b corresponding to the following constraints:
% sum(t(n)(j)(l),l=1:J)=1, any n, j
A1 = sparse(N*J, dim_x); b1 = ones(N*J,1);
for n = 1:N
    for j = 1:J
        row = (n-1)*J + j;
        for l = 1:J
            idx_exp_cone = ((n-1)*J+j-1)*J+l;
            idx_t_njl = (J+Q)+ N + 1 + (3*idx_exp_cone-1);
            A1(row, idx_t_njl) = 1;
        end
    end
end

% Assemble part of the matrix A in Ax=b corresponding to the following constraints:
% u(n)(j)(l) = 1, any n, j, l
A2 = sparse(N*J^2, dim_x); b2 = ones(N*J^2,1);
for n = 1:N
    for j = 1:J
        for l = 1:J
            row = ((n-1)*J+j-1)*J+l;
            idx_exp_cone = ((n-1)*J+j-1)*J+l;
            idx_u_njl = (J+Q)+ N+1 + 3*idx_exp_cone;
            A2(row, idx_u_njl) = 1;
        end
    end
end

% Assemble part of the matrix A in Ax=b corresponding to the s_njl constraints
A3 = sparse(N*J^2, dim_x); b3 = zeros(N*J^2, 1);
for n = 1:N
    for j = 1:J
        for l = 1:J
            row = ((n-1)*J+j-1)*J+l;
            idx_exp_cone = ((n-1)*J+j-1)*J+l;
            idx_s_njl = (J+Q)+ N + 1 + (3*idx_exp_cone-2);
            A3(row, idx_s_njl) = -1;
            idx_a_n = (J+Q)+n;
            A3(row, idx_a_n) = 1;
            idx_b = (J+Q)+N+1;
            A3(row, idx_b) = (j==j_hat(n));
            if l~=j
                idx_alpha_l = l; 
                idx_alpha_j = j;
                idx_beta = J+(1:Q);
                A3(row, idx_alpha_l) = 1;
                A3(row, idx_alpha_j) = -1;
                x_nl = [disp_mat(n,l); feat_mat(n,l); price_mat(n,l)];
                x_nj = [disp_mat(n,j); feat_mat(n,j); price_mat(n,j)];
                A3(row, idx_beta) = (x_nl - x_nj)';
            end
        end
    end
end

% alpha_var(1) = 0;
A4 = sparse(1, dim_x); b4 = 0;
idx_alpha_1 = 1;
A4(1, idx_alpha_1) = 1;

Abig = [A1; A2; A3; A4]; b = [b1; b2; b3; b4];

% Construct the cell array A
A{1} = Abig(:, 1:J);
A{2} = Abig(:, J+(1:Q));
A{3} = Abig(:, (J+Q)+(1:N));
A{4} = Abig(:, (J+Q)+N+1);
A{5} = Abig(:, (J+Q)+N+1+1:end);

disp('Done constructing A, c and b. Now calling the solver...');

% Fix training and testing datasets, vary Gamma and compare goodness of fit
clear Gamma;
% Initialize output arrays
opt_obj_max = zeros(length(Gamma_list),1);
likelihood_value = zeros(length(Gamma_list),1);
for idx = 1:length(Gamma_list)
    % Vary Gamma, which only appears in the objective function
    Gamma = Gamma_list(idx);
    c{4} = -(N - Gamma);
    input_options = [];
    if idx > 1
        clear input_options;
        input_options.initial_x = xre;
    end
    [opt_val, xre, yre, zre, info] = hsd_lqeu_fast(blk, A, c, b, input_options);
    % Get optimal alpha, beta, a[n], b and store them
    alpha_opt{idx} = xre{1};
    beta_opt{idx} = xre{2};
    a_opt{idx} = xre{3};
    b_opt{idx} = xre{4};
    opt_obj_max(idx) = -opt_val(1);
    % Perform validation
    disp('Begin validation...');
    disp('Given the estimates alpha_opt and beta_opt, compute the likelihood of the testing dataset...');
    % Now construct the N-by-J matrix representing this customer's purchase reord
    Ytest = purchase_record_mat_test(:,end-J+1:end);
    Ntest = size(purchase_record_mat_test,1);
    disp_mat_test = purchase_record_mat_test(:,2:1+J);
    feat_mat_test = purchase_record_mat_test(:,2+J:1+2*J);
    price_mat_test = purchase_record_mat_test(:,2+2*J:1+3*J);
    likelihood_value(idx) = 0;
    for n = 1:Ntest
        for j = 1:J
            x_n_1toJ = [disp_mat_test(n,:); feat_mat_test(n,:); price_mat_test(n,:)];
            x_nj = [disp_mat_test(n,j); feat_mat_test(n,j); price_mat_test(n,j)];
            curr_alpha = alpha_opt{idx};
            curr_beta = beta_opt{idx};
            expsum = sum(exp(curr_alpha-curr_alpha(j)+(x_n_1toJ-repmat(x_nj,[1,4]))'*curr_beta));
            likelihood_value(idx) = likelihood_value(idx) - Ytest(n,j)*log(expsum);
        end
    end
end
% Plot opt_obj_max and likelihood_value vs. Gamma_ratio
% plot(Gamma_ratio_list, opt_obj_max); hold on; plot(Gamma_ratio_list, likelihood_value); hold off;
plot(Gamma_list, likelihood_value);
title('Varying \Gamma');
xlabel('\Gamma'); 
%legend('optimal objective value (computed from the training dataset)', 'likelihood value (computed from the testing dataset)'); 
legend('likelihood value (computed from the testing dataset)'); 
legend('Location','southoutside');
legend('show');