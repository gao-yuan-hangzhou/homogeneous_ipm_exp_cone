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
spdata_mat = zeros(6*Ne,3); % Later, do spconvert(H_half_data_mat)

% Let input_vec = [u1;v1;w1; u2;v2;w2; ...]
all_u = input_vec(3*(1:Ne)-2);
all_v = input_vec(3*(1:Ne)-1);
all_w = input_vec(3*(1:Ne));
all_r = log(-all_v./all_u);
all_t = all_u - all_w + all_u.*all_r;

all_M11 = 1./(all_u.*all_t) + (all_r.^2)./(all_t.^2) + 1./(all_u.^2);
all_M12 = all_u.*all_r./(all_v.*(all_t.^2)) - 1./(all_v.*all_t); 
all_M13 = -all_r./(all_t.^2);

all_M22 = 1./(all_v.^2) + all_u./((all_v.^2).*all_t) + all_u.^2./((all_v.^2).*(all_t.^2)); 
all_M23 = -all_u./(all_v.*(all_t.^2));
all_M33 = 1./(all_t.^2);

% Form the columns of spdata_mat
col1 = zeros(6*Ne,1); col2 = zeros(6*Ne,1); col3 = zeros(6*Ne,1);
% Column 1 is [1;1;1;2;2;3; 4;4;4;5;5;6; ...]
col1(6*(1:Ne)-5) = (3*(1:Ne)-2)'; 
col1(6*(1:Ne)-4) = (3*(1:Ne)-2)'; 
col1(6*(1:Ne)-3) = (3*(1:Ne)-2)';
col1(6*(1:Ne)-2) = (3*(1:Ne)-1)';
col1(6*(1:Ne)-1) = (3*(1:Ne)-1)';
col1(6*(1:Ne)) = (3*(1:Ne))';
% Column 2 is [1;2;3;2;3;3; 4;5;6;5;6;6; ...]
col2(6*(1:Ne)-5) = (3*(1:Ne)-2)';
col2(6*(1:Ne)-4) = (3*(1:Ne)-1)';
col2(6*(1:Ne)-3) = (3*(1:Ne))';
col2(6*(1:Ne)-2) = (3*(1:Ne)-1)';
col2(6*(1:Ne)-1) = (3*(1:Ne))';
col2(6*(1:Ne)) = (3*(1:Ne))';
% Column 3 contain the actual values for the nonzero entries
col3(6*(1:Ne)-5) = all_M11/2;
col3(6*(1:Ne)-4) = all_M12;
col3(6*(1:Ne)-3) = all_M13;
col3(6*(1:Ne)-2) = all_M22/2;
col3(6*(1:Ne)-1) = all_M23;
col3(6*(1:Ne)) = all_M33/2;

spdata_mat = [col1, col2, col3];


% for k=1:Ne
%     u = input_vec(3*k-2); v = input_vec(3*k-1); w = input_vec(3*k); 
%     r = log(-v/u); t = u-w+u*r;
%     % Evaluate the Hessian of f_tilde_dual,
%     % The (second) barrier function for the dual of the exponential cone
%     M11 = 1/(u*t) + r^2/t^2 + 1/u^2; M12 =  u*r/(v*t^2) - 1/(v*t); M13 = -r/t^2;
%     M22 = 1/v^2 + u/(v^2*t) + u^2/(v^2*t^2); M23 = -u/(v*t^2);
%     M33 = 1/t^2;
%     spdata_mat(6*k-5:6*k,:) = [3*k-2, 3*k-2, M11/2; 3*k-2, 3*k-1, M12; 3*k-2, 3*k, M13;
%                                3*k-1, 3*k-1, M22/2; 3*k-1, 3*k, M23; 
%                                  3*k,   3*k, M33/2];
% end

H_half = spconvert(spdata_mat);
H = H_half+H_half';
end


