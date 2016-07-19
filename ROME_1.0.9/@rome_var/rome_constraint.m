function constraint_index = rome_constraint(obj)

% ROME_VAR\ROME_CONSTRAINT Adds a constraint into the current model
%
%   constraint_index = rome_constraint(obj)
%
%   obj must be a rome_var
%
% Modification History: 
% 1. Joel 

% check that the object is valid and remove degenerate case
if(isempty(obj) || obj.Cone == rome_constants.NO_CONE || obj.IsConst) 
    return;
end

% Get handle to the current model
h_curr_model = rome_get_current_model();

% Vectorize object here
S.type = '()';
S.subs = {':'};
if(obj.Cone ~= rome_constants.SOC)
    obj = obj.subsref(S);
end

% define output
constraint_index = [];

if(obj.IsCertain())
    % DETERMINISTIC VARIABLES
    % -----------------------
    constraint_index = rome_constraint_index;
    
    % represents Equality with ZERO
    if(obj.Cone == rome_constants.ZERO) 
        % append to linear constraints
        num_LCs = size(h_curr_model.LC, 1);
        h_curr_model.IndEq = [h_curr_model.IndEq; num_LCs + (1:obj.TotalSize)'];
        h_curr_model.LC    = [h_curr_model.LC; obj];

        % define constraint index output
        constraint_index.Index = (h_curr_model.LC.TotalSize - obj.TotalSize) + (1:obj.TotalSize).';
        constraint_index.Type = rome_constraint_index.LC;

    % represents Nonnegative Orthant Cone
    elseif(obj.Cone == rome_constants.NNOC) 
        % if non-diagonal
        if(~obj.IsDiag)
            h_curr_model.LC    = [h_curr_model.LC; obj];    % employ overloaded vertcat

            % define constraint index output
            constraint_index.Index = (h_curr_model.LC.TotalSize - obj.TotalSize) + (1:obj.TotalSize).';
            constraint_index.Type = rome_constraint_index.LC;

        % if diagonal, and positive multiplier
        elseif(obj.DiagMult > 0)
            h_curr_model.add_lower_bound(obj.NumUnmappedVars, -obj.BiAffineMap(:, 1) ./ obj.DiagMult );

            % define constraint index output
            constraint_index.Index = (obj.NumUnmappedVars) + (1:obj.TotalSize).';
            constraint_index.Type = rome_constraint_index.LB;

        % if diagonal, and negative multiplier
        elseif(obj.DiagMult < 0)
            h_curr_model.add_upper_bound(obj.NumUnmappedVars, -obj.BiAffineMap(:, 1) ./ obj.DiagMult);

            % define constraint index output
            constraint_index.Index = (obj.NumUnmappedVars) + (1:obj.TotalSize).';
            constraint_index.Type = rome_constraint_index.UB;
        else
            error('rome_var:rome_constraint:ZeroDiagMult', ...
                'Diagonal Multiplier is zero. Your rome_var is actually a constant and should not be put in a constraint!');
        end

    % represents Second Order Cone
    elseif(obj.Cone == rome_constants.SOC)
        % temp
%         if(1)
% Check if we have SQ non-linearity
        if(strcmp(obj.NonLinearity, 'SQ'))
            h_curr_model.QC = [h_curr_model.QC; {obj}];
        else
            % find the current number of variables
            old_num_vars = h_curr_model.NumVars;

            % find size of x
            sz_x = size(obj);
            nskip = sz_x(1);
            num_soc = prod(sz_x(2:end));  % number of embedded SOC constraints

            % For SOC, delay vectorization to here (S is defined at top)
            obj = obj.subsref(S);

            % construct the slack variables with their bounds
            x_s = rome_model_var(length(obj));
            h_curr_model.add_lower_bound(x_s.NumUnmappedVars, zeros(num_soc, 1), 0, nskip);

            % subtract the slack
            y = obj - x_s;

            % append to linear equality constraints
            num_LCs = size(h_curr_model.LC, 1);
            h_curr_model.IndEq = [h_curr_model.IndEq; num_LCs + (1:y.TotalSize)'];
            h_curr_model.LC    = [h_curr_model.LC; y];

            % define constraint index output (NOT YET)
            %         constraint_index.Index = (h_curr_model.LC.TotalSize - obj.TotalSize) + (1:obj.TotalSize).';
            %         constraint_index.Type = rome_constraint_index.LC;
            new_QC = 1:nskip:length(obj);
            new_QC = old_num_vars + [new_QC; new_QC + (nskip - 1)]';
            new_QC = mat2cell(new_QC, ones(1, size(new_QC, 1)), 2);
            h_curr_model.QC = cat(1, h_curr_model.QC, new_QC);        
        end
    end
    
elseif (obj.IsRand())
    % UNCERTAIN VARIABLES
    % -----------------    
    % represents Equality with ZERO
    if(obj.Cone == rome_constants.ZERO)
        % append to linear constraints
        num_LCs = size(h_curr_model.rndLC, 1);
        h_curr_model.rndIndEq = [h_curr_model.rndIndEq; num_LCs + (1:obj.TotalSize)'];
        h_curr_model.rndLC    = [h_curr_model.rndLC; obj];

        % if the object is a diagonal one
        if(obj.IsDiag)
            insert_vals = [obj.NumUnmappedRandVars + (1:obj.NumMappedRandVars)', -obj.BiAffineMap(:, 1) ./ obj.DiagMult];
            h_curr_model.rndEq = [h_curr_model.rndEq; insert_vals];
        else
%             error('rome_var:rome_constraint:RandEqProjection', 'Not implemented yet');
        end
        % represents Nonnegative Orthant Cone
    elseif(obj.Cone == rome_constants.NNOC)
        % if non-diagonal
        if(~obj.IsDiag)
            h_curr_model.rndLC  = [h_curr_model.rndLC; obj];    % employ overloaded vertcat

        % if diagonal, and positive multiplier
        elseif(obj.DiagMult > 0)
            h_curr_model.add_lower_bound(obj.NumUnmappedRandVars, -obj.BiAffineMap(:, 1) ./ obj.DiagMult, 1);

        % if diagonal, and negative multiplier
        elseif(obj.DiagMult < 0)
            h_curr_model.add_upper_bound(obj.NumUnmappedRandVars, -obj.BiAffineMap(:, 1) ./ obj.DiagMult, 1);
        else
            error('rome_var:rome_constraint:ZeroDiagMult', ...
                'Diagonal Multiplier is zero. Your rome_var is actually a constant and should not be put in a constraint!');
        end

        % represents Second Order Cone
    elseif(obj.Cone == rome_constants.SOC)
        %         error('rome_rand_var:rome_constraint:SOCNotAllowed', 'SOC
        %         Constraint not allowed for uncertain variables.');
        h_curr_model.rndQC = [h_curr_model.rndQC; {obj}];
    end
else
    % LDR Variables
    % -------------    
    % number of original linear inequality constraints prior to adding
    num_LCs = size(h_curr_model.ldrLC, 1);
    
    % Check for mean constraints
    start_ind = obj.NumUnmappedRandVars + 1;
    end_ind   = start_ind + obj.NumMappedRandVars - 1;
    test_vec = h_curr_model.ZIsMean(start_ind:end_ind);
    
    % is_mean is true for mean constraints
    if(~any(test_vec))
        is_mean = false; 
    else
        test_mat = repmat([false; test_vec]', obj.NumMappedVars + 1, 1);
        test_mat = test_mat(:)';
        is_mean = ~isempty(find(obj.BiAffineMap(:, test_mat), 1));
    end

    % represents equality with zero    
    if(obj.Cone == rome_constants.ZERO)
        h_curr_model.ldrIndEq = [h_curr_model.ldrIndEq; num_LCs + (1:obj.TotalSize)'];
        h_curr_model.ldrLC = [h_curr_model.ldrLC; obj];
    % represents Nonneg Orthant cone
    elseif(obj.Cone == rome_constants.NNOC)        
        if(is_mean)
            h_curr_model.ldrMeanInd = [h_curr_model.ldrMeanInd; num_LCs + (1:obj.TotalSize)'];
        end        
        h_curr_model.ldrLC = [h_curr_model.ldrLC; obj];
    % represents Second Order Cone
    elseif(obj.Cone == rome_constants.SOC)
        error('rome_var:rome_constraint:NotYetImplement', 'This constraint type has not been implemented for LDR yet');
    end

end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
