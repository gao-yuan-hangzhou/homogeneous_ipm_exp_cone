function f_val = f_exp_tilde(input_vec)
% input_vec is in (K_exp)^N, which is of dimension 3*N
N = length(input_vec)/3;
f_val = 0;
for k=1:N
    x = input_vec(3*k-2); y = input_vec(3*k-1); z = input_vec(3*k);
    omega_bar = wrightOmega(1-x/z-log(z)-log(y));
    f_val = f_val + (-2*log(z) - log(y)-log((1-omega_bar)^2/omega_bar) - 3);
end
end

