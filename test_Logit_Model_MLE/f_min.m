function f_val = f_min(a, b)
% Input: alpha (J-dimensional), beta (Q-dimensional)
% Output: objective value
% f_min is the negaive of the likelihood function

global Y N J Q disp_mat feat_mat price_mat;

g = @(alpha, beta, n, j) ...
    log(...
    sum(...
    exp(...
    (alpha-alpha(j))+([disp_mat(n,:)', feat_mat(n,:)', price_mat(n,:)']-repmat([disp_mat(n,j), feat_mat(n,j), price_mat(n,j)],[J,1]))*beta)));

% Sum up all terms y_nj * p_nj(alpha,beta)
f_val = 0;
for n = 1:N
    for j = 1:J
        f_val = f_val + Y(n,j)*g(a,b,n,j);
    end
end

% Function ends
end