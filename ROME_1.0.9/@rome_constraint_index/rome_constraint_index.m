% ROME_CONSTRAINT_INDEX Indexing system into a constraint in ROME
%
% Modification History: 
% 1. Joel 
% 

classdef rome_constraint_index
    % CONSTANTS
    % ----------
    properties (Constant = true, Hidden = true)
        % Constants representing Constraint Type
        LB  = -1;  % Lower Bound Constraint x >= L
        LC  =  0;  % Linear Constraint Ax (<=, ==, >=) b
        UB  =  1;  % Upper Bound Constraint x <= U        
    end

    properties 
        % Constraint Type: LB / LC / UB
        Type  = rome_constraint_index.LC; 
        
        % Actual Indexing into rows of model constraints
        Index = [];
    end
    
    methods
        % constructor
        function obj = rome_constraint_index()
        end
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
