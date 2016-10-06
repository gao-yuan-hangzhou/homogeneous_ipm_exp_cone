% load('likelihood_value.mat');
load('data_and_results_3291_samlpe.mat');

% Compute likelihood measure for the testing sets
% through both Method 1 and Method 2
num_split = length(alpha_opt);
Ntotal = size(purchase_record_mat,1);
J = 4;
disp('Begin validation...');
disp('Given the estimates alpha_opt and beta_opt, compute the likelihood of the testing dataset...');
for dataset_index=1:num_split
    purchase_record_mat_test = purchase_record_mat(floor(Ntotal/num_split)*(dataset_index-1)+1:floor(Ntotal/num_split)*dataset_index,:);
    % Now construct the N-by-J matrix representing this customer's purchase reord
    Y = purchase_record_mat_test(:,end-J+1:end);
    N = size(purchase_record_mat_test,1);
    disp_mat = purchase_record_mat_test(:,2:1+J);
    feat_mat = purchase_record_mat_test(:,2+J:1+2*J);
    price_mat = purchase_record_mat_test(:,2+2*J:1+3*J);
    likelihood_value{dataset_index} = 0;
    meam_square_residual_value{dataset_index} = 0;
    for n=1:N
        for j=1:J
            x_n_1toJ = [disp_mat(n,:); feat_mat(n,:); price_mat(n,:)];
            x_nj = [disp_mat(n,j); feat_mat(n,j); price_mat(n,j)];
            curr_alpha = alpha_opt{dataset_index};
            curr_beta = beta_opt{dataset_index};
            expsum = sum(exp(curr_alpha-curr_alpha(j)+(x_n_1toJ-repmat(x_nj,[1,4]))'*curr_beta));
            likelihood_value{dataset_index} = likelihood_value{dataset_index} - Y(n,j)*log(expsum);
            p_hat_n_j = expsum;
            meam_square_residual_value{dataset_index} = meam_square_residual_value{dataset_index} + (1/N)*(Y(n,j)-p_hat_n_j)^2;
        end
    end
end

opt_obj_max
likelihood_value
meam_square_residual_value