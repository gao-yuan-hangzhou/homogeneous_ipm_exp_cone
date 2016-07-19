function varargout = display(sol_obj, bDebug)
%
% ROME_SOL\DISPLAY Prettyprints object printing it nicely
% in a string. This assumes the object is a vector. If the object is a
% matrix or N-dimensional, you will have to reshape the output on your own.
%
% USAGE:
% str = x.display;  % outputs prettyprint into a string
% x                 % will display in the MATLAB command window
% 
%
% FOR DEVELOPERS:
% display(x, true) % can see underlying structure of x
% 
% History
% 1. Created by Joel 14 May 2009 
%

% HACK For developers
if(nargin == 2)
    builtin('disp', sol_obj);
    return;
end

% instantiate stringc
str = sprintf('\n%s = \n\n', inputname(1));

% extract values
x_val = sol_obj.LDRAffineMap;
x_coeff = sol_obj.DeflectCoeff;
x_deflect_vals = sol_obj.DeflectAffineMap;

% build strings
sign_str = '-+';
abs_x_val = abs(x_val);
sign_x_val = sign_str(1 + ceil((sign(x_val) + 1)/2));

abs_x_coeff = abs(x_coeff);
sign_x_coeff = sign_str(1 + ceil((sign(x_coeff) + 1)/2));


abs_x_deflect_vals = abs(x_deflect_vals);
sign_x_deflect_vals  = sign_str(1 + ceil((sign(x_deflect_vals) + 1)/2));

% for each output variable
for ii = 1:size(x_val, 1)
    % for each component of primitive uncertainty
    for jj = 1:size(x_val, 2)
        if(jj == 1)
            str = [str, sprintf('% 0.3f', x_val(ii, jj))];
        else
            if(x_val(ii,jj) == 0)
                continue;
            end
            str = [str, sprintf(' %s %0.3f*z%d', sign_x_val(ii, jj), abs_x_val(ii, jj), jj - 1)];
        end
    end
    
    % append nonlinear deflected components
    for jj = 1:size(x_coeff, 2)
        if(x_coeff(ii, jj) ~= 0)
            str = [str, sprintf(' %s %0.2f(', sign_x_coeff(ii, jj), abs_x_coeff(ii, jj))];
            for kk = 1:size(x_deflect_vals, 2)
                if(x_deflect_vals(jj,kk) == 0)
                    continue;
                end
                if(kk == 1)
                    str = [str, sprintf('%0.3f', x_deflect_vals(jj, kk))];
                else
                    str = [str, sprintf(' %s %0.3f*z%d', sign_x_deflect_vals(jj, kk), abs_x_deflect_vals(jj, kk), kk - 1)];
                end                
            end
            str = [str, ')^-']; 
        end
    end
    str = [str, sprintf('\n')];
end

if(nargout == 0)
    disp(str);
else
    varargout(1) = str;
end



% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
