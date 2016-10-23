function output = is_in_exp_cone_interior(x)
% Check whether a given vector is in (K_exp)^N
N = length(x)/3;
if floor(N) ~= N
    display('Warning: dimension of x is not a multiple of 3!');
else
    output = true;
    for k=1:N
        x3dim = x(3*k-2:3*k);
        if ~(x3dim(3)>0 && x3dim(2)>0 && exp(x3dim(1)/x3dim(3)) < x3dim(2)/x3dim(3))
            output = false;
            break;
        end
    end
end