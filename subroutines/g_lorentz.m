function g = g_lorentz(x)
% Compute the gradient of the barrier for the Lorentz cone at x in Q(n),
% where Q(n) = {(x(1),...,x(n))| norm(x(2:n))<=x(1)} 
% We choose the log barrier with v = 1 as given in CVX Guide: http://cvxr.com/cvx/doc/CVX.pdf
g = [-x(1); x(2:end)]/(x(1)^2 - x(2:end)'*x(2:end));
end

