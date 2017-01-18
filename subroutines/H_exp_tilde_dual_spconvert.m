function H = H_exp_tilde_dual_spconvert(input_vec)
% Compute the Hessian of the dual exp barrier (from the second pair)
% where N is the number of dual exp cones
Ne = length(input_vec)/3;
if floor(Ne)~=Ne
    display('Error: in g_exp_tilde_dual, dim(input_vec) is not a multiple of 3!');
end
H = sparse(3*Ne,3*Ne);
% There are 6*Ne distinct entries in H_exp
% We find HH = spconvert(H_data_mat)
% such that H_exp = (HH+HH')/2
H_data_mat = zeros(6*Ne,3); % Later, do spconvert(H_data_mat)
for k=1:Ne
    u = input_vec(3*k-2); v = input_vec(3*k-1); w = input_vec(3*k); 
    r = log(-v/u); t = u-w+u*r;
    % Evaluate the Hessian of f_tilde_dual,
    % The (second) barrier function for the dual of the exponential cone
    M11 = 1/(u*t) + r^2/t^2 + 1/u^2; M12 =  u*r/(v*t^2) - 1/(v*t); M13 = -r/t^2;
    M22 = 1/v^2 + u/(v^2*t) + u^2/(v^2*t^2); M23 = -u/(v*t^2);
    M33 = 1/t^2;
    H_data_mat(6*k-5:6*k,:) = [3*k-2, 3*k-2, M11/2; 3*k-2, 3*k-1, M12; 3*k-2, 3*k, M13;
                               3*k-1, 3*k-1, M22/2; 3*k-1, 3*k, M23; 
                                 3*k,   3*k, M33/2];
end

H_half = spconvert(H_data_mat);
H = H_half+H_half';
end


