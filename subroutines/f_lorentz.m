function f_val = f_lorentz(x)
% Evaluate the barrier for the Lorentz cone
% Q(n) = {(x(1),...,x(n))| norm(x(2:n))<=x(1)} 
f_val = -(1/2) * (log(x(1)^2 - x(2:end)'*x(2:end)));
end