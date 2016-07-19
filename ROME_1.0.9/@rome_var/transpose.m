function out_var_obj = transpose(A)

% ROME_VAR\TRANSPOSE Implements complex conjugate transpose in ROME
%
%   C = A.'
%
%   Same as complex transposition
%
% Modification History: 
% 1. Joel 

% Throw an error if not a matrix
if(length(size(A)) > 2)
    error('rome_var:transpose:InvalidArg', 'Transpose only operates on 2 dimensional matrices');
end

out_var_obj = A';


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
