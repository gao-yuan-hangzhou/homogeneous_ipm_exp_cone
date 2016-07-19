function out_var_obj = eq(A, B)

% ROME_VAR\LE Implements equal to operator on rome_var objects
%
%   A == B 
%
%   At least one of A, B must be a rome_constants.
%
%   Returns a new rome_var object, B-A, constrained to be zero. Also
%   registers the constraint with the current model
%
% Modification History: 
% 1. Joel 

% Not necessary to error check for size matching, since we will handle it
% implicitly upon calling the 'minus' operator

out_var_obj = B - A;
out_var_obj.Cone = rome_constants.ZERO;

% % registers constraint with model
% S.type = '()';
% S.subs = {':'};
% rome_constraint(out_var_obj.subsref(S));


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
