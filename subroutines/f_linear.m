function f_val = f_linear(x)
% Evaluate the barrier for the nonnegative orthant
f_val = -sum(log(x));
end