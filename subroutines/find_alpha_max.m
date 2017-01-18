function alpha_max = find_alpha_max(x_bar, predictor_search_direction, dimension_info)
% x_bar = [x, y, z, tau, kappa, theta]
% Finds the maximum step length alpha_max such that
% x + alpha_max * dx is in int(K),
% z + alpha_max * dz is in int(K^*),
% tau + alpha_max * dtau > 0,
% kappa + alpha_max * dkappa > 0.
% If the above inequalities are satisfied? one has theta + alpha_max * dtheta > 0

% Extract x, y, z, tau, kappa
% x = x_bar(1:Nt); y = x_bar(Nt+1:Nt+m); z = x_bar(Nt+m+1:2*Nt+m);
% tau = x_bar(2*Nt+m+1); kappa = x_bar(2*Nt+m+2); theta = x_bar(2*Nt+m+3);
% Extract dx, dy, d, dtau, dkappa
% dx = predictor_search_direction(1:Nt); dy = predictor_search_direction(Nt+1:Nt+m); dz = predictor_search_direction(Nt+m+1: 2*Nt+m); 
% dtau = predictor_search_direction(2*Nt+m+1); dkappa = predictor_search_direction(2*Nt+m+2); dtheta = predictor_search_direction(2*Nt+m+3);

% Bisection search
alpha_low=0; alpha_high = 1; max_it = 20;
for it_count = 1:max_it
    alpha_middle = 0.5*alpha_low + 0.5*alpha_high;
    if is_in_product_cone(x_bar + alpha_high*predictor_search_direction, dimension_info)
        alpha_low = alpha_middle;
    else
        alpha_high = alpha_middle;
    end
    if (it_count >= 5 && alpha_low >= 0.1)
        break
    end
end
alpha_max = alpha_low;
end
