function v = v_dual(z, dimension_info)
% Compute the vector v such that H_dual(z) = E_dual + v*v'
if length(z) ~= dimension_info.l+sum(dimension_info.q) + 3*dimension_info.e
    display('Error: In H_dual(z, dimension_info), length of x does NOT match dimension_info!');
end

Nl = dimension_info.l; 
Nq = dimension_info.q; 
Ne = dimension_info.e;

vq = zeros(sum(Nq), 1);
for k = 1:length(Nq)
    vq(1+sum(Nq(1:k-1)):sum(Nq(1:k))) = v_lorentz(z(Nl+1+sum(Nq(1:k-1)):Nl+sum(Nq(1:k))));
end

v = zeros(Nl+sum(Nq)+3*Ne, 1);
v(Nl+1:Nl+sum(Nq)) = vq;
end

