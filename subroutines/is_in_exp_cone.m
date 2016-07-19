function output = is_in_exp_cone(x)
% Check whether a given vector is in (K_exp)^N
tol = 0;
N = length(x)/3;
if floor(N) ~= N
    display('Warning: dimension of x is not a multiple of 3!');
else
    output = true;
    for k=1:N
        x_3dim = x(3*k-2:3*k);
        if ~((abs(x_3dim(3)) <= tol && x_3dim(1) <=tol && x_3dim(2) >= -tol) || exp(x_3dim(1)/x_3dim(3)) <= x_3dim(2)/x_3dim(3) + tol) 
            % display(['Check k=' num2str(k)]);
            % display(x_3dim);
            output = false;
            break;
        end
    end
end

