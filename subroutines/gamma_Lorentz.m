function output = gamma_Lorentz(z)
% Compute gamma(z) for z in Q(p)
% gamma(z) = sqrt(z(1) - z(2:end)'*z(2:end))
output = sqrt(z(1)^2 - z(2:end)'*z(2:end));
end

