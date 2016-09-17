addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set sample size, number of subsets and Gamma/N, where N is the size of the training set
sample_size = 1000; % sample_size = 100;
num_split = 5;
Gamma_ratio = 0.1;
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
% load the same purchase_record_mat as before so no resampling
% load('./cross_valudation_200_observations_5_subsets_Gamma_40/data_and_results_1000_samlpe.mat');

% N is the number of observations of this customer
Ntotal = size(purchase_record_mat,1);
disp(['Total number of observations = ' num2str(Ntotal)]);

% Q is the number of attributes
Q = 3;

dataset_length = floor(Ntotal/num_split);
for dataset_index=1:num_split
    % Partition the data into training and testing data
    purchase_record_mat_train = purchase_record_mat([1:floor(Ntotal/num_split)*(dataset_index-1), floor(Ntotal/num_split)*dataset_index+1:end],:);
    %purchase_record_mat_test = purchase_record_mat(floor(N/num_split)*(dataset_index-1)+1:floor(N/num_split)*dataset_index,:);
    % Now construct the N-by-J matrix representing this customer's purchase reord
    Y = purchase_record_mat_train(:,end-J+1:end);
    N = size(purchase_record_mat_train,1); Gamma = 0.01*N;
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
    c{4} = -(N - Gamma);
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

    % Call hsd_lqeu() or hsd_lqeu_Schur_dense_column_handling()
    % [opt_val, xre, yre, zre, info] = hsd_lqeu(blk, A, c, b);
    % [opt_val, xre, yre, zre, info] = hsd_lqeu_Schur_dense_column(blk, A, c, b);
    [opt_val, xre, yre, zre, info] = hsd_lqeu_fast(blk, A, c, b);
    %[opt_val, xre, yre, zre, info] = hsd_lqeu(blk, A, c, b);

    % Get optimal alpha, beta, a[n], b and store them
    alpha_opt{dataset_index} = xre{1};
    beta_opt{dataset_index} = xre{2};
    a_opt{dataset_index} = xre{3};
    b_opt{dataset_index} = xre{4};
    opt_obj_max{dataset_index} = -(opt_val(1)+opt_val(2))/2;
end

% Save all results
mkdir('cross_validation_data_and_result');
save('./cross_validation_data_and_result/data_and_results_1000_samlpe.mat', 'purchase_record_mat', 'alpha_opt', 'beta_opt', 'a_opt', 'b_opt', 'opt_obj_max');
% Perform validation
disp('Begin validation...');
disp('Given the estimates alpha_opt and beta_opt, compute the likelihood of the testing dataset...');
for dataset_index=1:num_split
    purchase_record_mat_test = purchase_record_mat(floor(Ntotal/num_split)*(dataset_index-1)+1:floor(Ntotal/num_split)*dataset_index,:);
    % Now construct the N-by-J matrix representing this customer's purchase reord
    Y = purchase_record_mat_test(:,end-J+1:end);
    N = size(purchase_record_mat_test,1); % Gamma = 0.1*N;
    disp_mat = purchase_record_mat_test(:,2:1+J);
    feat_mat = purchase_record_mat_test(:,2+J:1+2*J);
    price_mat = purchase_record_mat_test(:,2+2*J:1+3*J);
    likelihood_value{dataset_index}=0;
    for n=1:N
        for j=1:J
            x_n_1toJ = [disp_mat(n,:); feat_mat(n,:); price_mat(n,:)];
            x_nj = [disp_mat(n,j); feat_mat(n,j); price_mat(n,j)];
            curr_alpha = alpha_opt{dataset_index};
            curr_beta = beta_opt{dataset_index};
            expsum = sum(exp(curr_alpha-curr_alpha(j)+(x_n_1toJ-repmat(x_nj,[1,4]))'*curr_beta));
            likelihood_value{dataset_index} = likelihood_value{dataset_index} - Y(n,j)*log(expsum);
        end
    end
end

disp(' ');
disp('Optimal objective (maximum likelihood) of training sets = ');
disp([opt_obj_max{1},opt_obj_max{2}, opt_obj_max{3}, opt_obj_max{4}, opt_obj_max{5}]);
disp('Likelihood values of testing sets = ');
disp([likelihood_value{1},likelihood_value{2}, likelihood_value{3}, likelihood_value{4}]);
save('./cross_validation_data_and_result/likelihood_value.mat', 'likelihood_value');