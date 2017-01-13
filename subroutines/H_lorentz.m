function H = H_lorentz(x) % x is a column vector
% Compute the Hessian of the log barrier of the Lorentz cone 
% Q(m) = {x|sqrt(x(2)^2+...+x(n)^2)<=x(1)}
% THe log barrier of Q(m) chosen is 
% f(x) = -(1/2)*log(x(1)^2 - (x(2)^2+...+x(n)^2)) with v = 1
xTJx = x(1)^2 - x(2:end)'*x(2:end);
H = xTJx^(-2) * (2*[x(1);-x(2:end)]*[x(1); -x(2:end)]' - xTJx*diag([1; -ones(length(x)-1,1)]));
end