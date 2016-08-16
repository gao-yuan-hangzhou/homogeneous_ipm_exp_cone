function z_inv = Lorentz_inv(z)
% Compute z^(-1) where z is in the lorentz cone Q(p)
z_inv = [z(1); -z(2:end)]/(z(1)^2 - z(2:end)'*z(2:end));
end

