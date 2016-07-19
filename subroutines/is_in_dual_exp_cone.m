function output = is_in_dual_ezp_cone(z)
tol = 0;
N = length(z)/3;
if floor(N) ~= N
    display('Warning: dimension of z is not a multiple of 3!');
else
    output = true;
    for k=1:N
        z_3dim = z(3*k-2:3*k);
        if ~((abs(z_3dim(1)) <= tol && z_3dim(2) >= -tol && z_3dim(3) >= -tol) || exp(z_3dim(3)/z_3dim(1)) <= -exp(1)*z_3dim(2)/z_3dim(1) + tol) 
            output = false;
            break;
        end
    end
end
