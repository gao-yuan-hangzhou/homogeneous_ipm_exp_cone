function [X_opt, primobj, iteration_count,res] ...
          = badmm_assignment_problem_solver(C,a,b, rho, total_number_of_iterations)
% Solves an assignment problem by taking in the cost matrix C and
% right-hand-side vectors a and b.
% C is m-by-n, a is m-dimensional, b is n-dimensional
% Inputs: C,a,b, rho, total_number_of_iterations
% Outputs: X_opt, optimal_obj_val, iteration_count,res_output, obj_val_array

% ADMM Formulation: min <C,X> s.t. X in S(X), Z in S(Z), X - Z = 0 (A = I, B = -I, c = 0)

[m,n]=size(C);
% Initialize X and Z: note that we must keep them positive

X = repmat(a/n, [1,n]);
Z = repmat(b'/m, [m,1]);

% The lagrange multiplier
Y = zeros(m,n);

% The main iteration of BADMM
res = zeros(total_number_of_iterations,1);
breakyes = 0; 
for k = 1:total_number_of_iterations
    % Updating X, row-major
    Temp = Z.*exp(-(C+Y)/rho);
%   x_row_sums = Temp * ones(n,1);
    X_new = spdiags(a./(Temp * ones(n,1) + eps),0,m,m)*Temp;
    
    % Updating Z, column-major
    Temp = X_new.*exp(Y/rho);
%   z_column_sums = ones(1,m)*Temp;
    Z_new = Temp * spdiags(b./(ones(1,m)*Temp  + eps)',0,n,n);
    
    if (true)
        % Calculate the residual, in this case gamma = 0.25
        res(k) = sum(sum(X_new.* log((X_new+eps)./(Z+eps)))) + 0.25*sum(sum((X_new - Z_new).^2));
        if (res(k) < 1e-5)
            breakyes = 1; 
        end
    end
             
        
    % Replace the old X and Z
    X = X_new;
    Z = Z_new;
    % Update Y, note taht we must take tau <= (1-2*gamma)*rho = 0.5*rho
    Y = Y + 0.5*rho * (X - Z);
        
    if (rem(k,50)==1) || (breakyes)
        primobj = sum(sum(C.*X_new)); 
        fprintf('\n k = %4.0f, res = %3.2e, %7.6e',k,res(k),primobj);
    end
    if (breakyes); break; end
end
X_opt = X;
iteration_count = k;

