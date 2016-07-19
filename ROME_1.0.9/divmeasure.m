function [x_mean, x_stddev, x_fdev,x_bdev,x_low,x_upp] = divmeasure(x, w, tol, M)
%DIVMEASURE Return deviation measures of x with weights w
%   [x_mean, x_stddev, x_fdev,bdev,x_low,x_upp] =  DIVMEASURE(x,w)
%       Inputs:
%           x data inputs 
%           w data weights: assume equal weights if w is absent
%           tol : tolerance, default = 1E-8
%           M   : large number, default = 100;
%
%       Outputs:
%           x_mean is the mean  
%           x_stddev is the standard deviation 
%           x_fdev is the forward deviation
%           x_bdev is the backward deviation
%           x_low is the lower support, max(x)
%           x_upp is the upper support, min(x)
%
% Last Modified:  
% 27 July 2009

% data length
x=x(:);
L = length(x);

% setting default arguments
if(nargin < 4)
    M = 100;
end
if(nargin < 3)
    tol = 1E-8;
end
if (nargin < 2)
    w = ones(L,1)/L;
end

% normalize
w = w(:)/sum(w);

%%%% Estimate Directional Deviations
x_mean = w'*x;
z = x - x_mean;
x_stddev = sqrt(sum(w.*z.^2));
x_upp = max(x);
x_low = min(x);

theta1 = tol;
theta4 = M/max(z);
theta2 =  theta1 + (theta4-theta1)/3;
theta3 =  theta1 + 2*(theta4-theta1)/3;

% Peform bisection search to find x_fdev
while theta4-theta1>tol
    tmp2 = sqrt(2*log(sum(exp(theta2*z).*w))/theta2^2);
    tmp3 = sqrt(2*log(sum(exp(theta3*z).*w))/theta3^2);
    %disp([theta2 tmp2 theta3 tmp3])
    if tmp2>tmp3
        theta4 = theta3;
    else
        theta1 = theta2;
    end;
    theta2 =  theta1 + (theta4-theta1)/3;
    theta3 =  theta1 + 2*(theta4-theta1)/3;
end;

x_fdev = max((tmp2+tmp3)/2,x_stddev);

theta1 = tol;
theta4 = -M/min(z);
theta2 =  theta1 + (theta4-theta1)/3;
theta3 =  theta1 + 2*(theta4-theta1)/3;

% Peform bisection search to find x_bdev,
while theta4-theta1>tol
    tmp2 = sqrt(2*log(sum(exp(-theta2*z).*w))/theta2^2);
    tmp3 = sqrt(2*log(sum(exp(-theta3*z).*w))/theta3^2);
    %disp([theta2 tmp2 theta3 tmp3])
    if tmp2>tmp3
        theta4 = theta3;
    else
        theta1 = theta2;
    end;
    theta2 =  theta1 + (theta4-theta1)/3;
    theta3 =  theta1 + 2*(theta4-theta1)/3;
end;
x_bdev = max((tmp2+tmp3)/2,x_stddev);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
