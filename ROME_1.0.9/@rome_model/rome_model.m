% ROME_MODEL
% DESCRIBES AN INSTANCE OF A MODEL IN ROME ACHITECTURE
%
%
% Modification History: 
% 1. Joel 
% 
classdef rome_model < handle

    properties (Constant = true, Hidden = true)
        MINIMIZE = -1;
        MAXIMIZE =  1;
    end
    
    % PRIVATE PROPERTIES
    % ------------------
    % Note to Joel: All properties are currently public
    properties % (SetAccess = protected, GetAccess = protected) % TO BE REMOVED LATER
        
        % General Properties
        Name  = [];      % Model Name
        MinMaxFlag = rome_model.MINIMIZE;    % Flag to set whether min or max (Default = minimize)
        ObjFn = [];      % Objective Function
        XSol  = [];      % Primal Solution
        ObjVal= NaN;     % Objective Function Value at Solution
        DualVars   = []; % Dual Variables
        DualObjVal = NaN;% Dual Objective Function
        ZRnd  = [];      % Uncertain Variables
        ZIsMean = [];    % Vector of flags indicating if the variable represents a mean 
        bBoundUsed = false; % scalar flag to signify that bound has been used
        
        % Properties of Deterministic Variables
        % -------------------------------------
        LC    = [];     % linear constraints
        QC    = {};     % Quadratic Constraints (TO BE CHANGED TO CONIC!)
        IndEq = [];     % Index of Equalities
        LB = [];        % Lower Bound Constraints
        UB = [];        % Upper Bound Constraints
        VarType = [];   % Variable Type ('C', 'I', 'B')
        
        % Properties of Uncertain Variables
        % -------------------------------------
        % 1. Support Information
        rndLC    = [];     % linear constraints
        rndQC    = {};     % Quadratic Constraints - NOT ALLOWED
        rndIndEq = [];     % Index of Equalities (linear)
        rndLB    = [];     % Lower Bound Constraints
        rndUB    = [];     % Upper Bound Constraints
        rndVarType = [];   % Variable Type ('C') - ONLY CONTINUOUS ALLOWED
        rndEq    = [];     % Index of PRIMITIVE variables which have equality constraints
        
        % 2. Parameter information
        rndMean  = [];     % Mean (vector)
        rndCovar = [];     % Covariance (vector or matrix)
        rndCovarMix = [];  % Mixing Matrix for Covariance (rand_var)
        rndFDev  = [];     % Forward Deviation (vector)
        rndBDev  = [];     % Backward Deviation (vector)
        rndDirDevMix = []; % Mixing Matrix for Directional Deviation (rand_var)
        
        % Properties of LDR Variables
        % -------------------------------------
        ldrVars  = [];     % LDR variables which have been declared (vectorized)
        ldrLC    = [];     % LDR constraints
        ldrIndEq = [];     % LDR Ind Equalities
        ldrMeanInd= [];    % Indices of Linear Expectation Constraints
        
        % Properties of Deflections
        % -------------------------------------
        ldrRemLC = [];     % LDRs removed because of deflections
        ldrDefCoeff = [];  % Matrix of Deflection coefficients
        ldrDefInd = [];    % Vector of deflection indices
        
        % Deflections Options
        % -------------------
        defObjFn       = [];    % Deflection Objective Function
        defLDRRestrict = [];    % Vector of restricted LDRs for deflection
        defBounds      = {};    % Cell array of fn handles to deflection bounds
        
        % Solver Choice
        Solver = '';
%       Solver='MOSEK';
%       Solver='SDPT3DUAL';
    end
    properties (Dependent = true)
        NumVars;
        NumRandVars;
    end
    
    % PROPERTY ACCESS METHODS
    % -----------------------
    methods
        function value = get.NumVars(obj)
            value = length(obj.XSol);
        end
        function value = get.NumRandVars(obj)
            value = length(obj.ZRnd);
        end

    end
    
    % METHODS
    % --------
    % Note: To be moved to separate files for better control 
    methods
        % Constructor
        % ------------
        function obj = rome_model(model_name)
            % global ROME Environment
            global ROME_ENV;
           
            if(nargin == 0 || isempty(model_name))
                % Create a dummy name if model is empty
                obj.Name = sprintf('rome_model_%0.4d', ROME_ENV.num_models);
            elseif(~ischar(model_name))
                % if model_name is not character array, error, and exit
                error('rome_model:ctor:InvalidArg', ...
                            'Constructor requires char array (string) argument');
            else
                obj.Name = model_name;                
            end
            
            % assign solver
            obj.Solver = ROME_ENV.DEF_SOLVER;
            
            % register this model with the environment
            rome_reg_model(obj);
        end
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
