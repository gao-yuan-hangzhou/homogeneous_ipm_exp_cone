function out_var_obj = max(obj, varargin)

% ROME_VAR\MAX Returns a rome_var object which is the maximum of two or more
% objects. Objects must be of the same size or scalars
%
% Usage:
%
%   C = max(A, B)
%   C = max(A, 2)
%
% Modification History: 
% 1. Joel 

% catch pathological case
if(nargin == 1)
    out_var_obj = obj;
    return;
end

% otherwise, we create a new object and enforce maximum
out_var_obj = rome_clone(obj);
rome_constraint(out_var_obj >= obj);

% enforce maximum on each guy
for ii = 1:length(varargin)
    rome_constraint(out_var_obj >= varargin{ii});
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
