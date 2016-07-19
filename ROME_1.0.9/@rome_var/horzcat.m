function out_var_obj = horzcat(varargin)

% ROME_VAR\HORZCAT Implements horizontal concatenation for rome_var
%
%   C = [A, B, ...]
%
%
% Modification History: 
% 1. Joel 

out_var_obj = cat(2, varargin{:});

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
