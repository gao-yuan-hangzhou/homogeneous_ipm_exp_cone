function output= is_in_lorentz_cone_interior(x)
% Check whether x has dimension at least 2 AND x is in the interior of a Lorentz cone
if length(x)<=1
    display('Error: x has dimension less than or equal to 1!');
    output = false;
else
    output = (x(1)>0) && (x(1)^2 > x(2:end)'*x(2:end));
end
end

