function f_val = f_exp_tilde_dual(input_vec)
% Evaluate the barrier for the dual exponential cone
% where input_vec is in (P_exp)^N, where N is the number of dual exp cones
N = length(input_vec)/3;
f_val = 0;
for k = 1:N
    u=input_vec(3*k-2); v=input_vec(3*k-1); w=input_vec(3*k);
    r=log(-v/u); t=u-w+u*r;
    f_val = f_val + (-log(-t) - log(-u) - log(v));
end

end

