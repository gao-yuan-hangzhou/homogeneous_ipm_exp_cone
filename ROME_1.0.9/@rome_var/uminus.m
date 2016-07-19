function A = uminus(A)

% ROME_VAR\UMINUS Implements unary minus in ROME
%
%   C = -A
%
%   A must be a rome_var object. Returns C, another rome_var object
%
% Modification History: 
% 1. Joel 

% Allocate space for output variable
A.BiAffineMap = -A.BiAffineMap;


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
