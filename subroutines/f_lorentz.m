function F = F_lorentz(x) % x is a column vector
% Compute W such that W*W = H_lorentz(x)
% Expression of H_sqrt can be found in
% Section 2.2 in http://www.seas.ucla.edu/~vandenbe/publications/coneprog.pdf
p = length(x);
xTJx = x(1)^2 - x(2:end)'*x(2:end);
M11 = x(1); 
M12 = -x(2:end)';
M21 = -x(2:end);
M22 = (x(2:end)*x(2:end)')/(x(1) + sqrt(xTJx)) + sqrt(xTJx)*eye(p-1);
F = 1/xTJx * [M11, M12; M21, M22]; 
end