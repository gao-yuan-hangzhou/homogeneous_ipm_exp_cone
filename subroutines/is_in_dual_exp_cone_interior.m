function output = is_in_dual_exp_cone_interior(z)
N = length(z)/3;
if floor(N) ~= N
    display('Warning: dimension of z is not a multiple of 3!');
else
    output = true;
    for k=1:N
        z3dim = z(3*k-2:3*k);
        if ~(z3dim(1)<0 && z3dim(2)>0 && exp(z3dim(3)/z3dim(1))<-exp(1)*z3dim(2)/z3dim(1))
            output = false;
            break;
        end
    end
end
