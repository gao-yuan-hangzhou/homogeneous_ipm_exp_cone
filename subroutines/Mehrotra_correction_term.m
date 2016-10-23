function T = Mehrotra_correction_term(pred_dir, dimension_info)
% Calculate the Mehrotra correction term for (dx_linear, dz_linear) and (dx_soc, dz_soc)

% Obtain dimension information
Nl = dimension_info.l; Nq = dimension_info.q; Ne = dimension_info.e; m = dimension_info.m; Nt = Nl + sum(Nq) + 3*Ne;
% Obtain dx_p and dz_p
dx_p = pred_dir(1:Nt); dz_p = pred_dir(Nt+m+1: 2*Nt+m);
% Set the return vector
T = zeros(Nt,1);
T(1:Nl) = dx_p(1:Nl) .* dz_p(1:Nl);    
end

