function output = precond_M_inv(h_hat, L_chol, U4, LG, UG, Msp_inv_U)
% The preconditioner function for the BiCGSTAB subroutine
% Given L = chol(Msp), [LG, UG] = lu(G), Msp_inv_U = Msp\U = L\(L'\U)
% Compute output = M\rhs approximately by
% output = L\(L'\h_hat) - Msp_inv_U*(UG\(LG\(U'*(L\(L'\h_hat)))))

Msp_inv_chol = @(X) L_chol\(L_chol'\X);
Msp_inv_h_hat = Msp_inv_chol(h_hat);
output = Msp_inv_h_hat - Msp_inv_U*(UG\(LG\(U4'*Msp_inv_h_hat)));
%disp('precond function is called!');

end

