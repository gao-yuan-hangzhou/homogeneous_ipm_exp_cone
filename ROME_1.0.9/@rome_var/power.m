function y = power(x, n)

% ROME_VAR\POWER Returns y such that x.^n <= y
%
% Currently only accepts scalar, positive integer n. x should be >= 0 in
% general.
%
% Modified By:
% 1. Joel

y = powercone(x, 1, n);

% if(~isscalar(n) || floor(n) ~= n || n < 1)
%     error('rome_var:power:NeedPosIntScalarPower', 'Power must be a positive integer');
% end
% 
% % degenerate case
% if(n == 1)
%     y = x;
%     return;
% end
% 
% % make output variable
% y = rome_model_var('Cone', rome_constants.NNOC);
% 
% % define starting points of iteration
% rhs = y;
% lhs = x;
% 
% while(n > 1)
%     if(mod(n,2) == 0)
%         % make new sq constraint 
%         lhs = sq(lhs);
%     else
%         % make new hypoquad constraint
%         rhs = hypoquad(lhs, rhs);        
%     end
%     
%     % update n
%     n = ceil (n / 2);
% end
% 
% % final 
% rome_constraint(lhs <= rhs);


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
