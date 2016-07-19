function [prob,add_obj] = convert_cplex(model_obj)

%
% ROME_MODEL\CONVERT_CPLEX Converts the Model object into cplex solvable
% form
% 

% Using CPLEX:
%
% x = cplexmiqcp(prob)
%
% Description:
%
% Finds the minimum of a problem specified by
% min      0.5*x'*H*x+f*x or f*x
% st.      Aineq*x      <= bineq
%          Aeq*x         = beq
%          l*x + x'*Q*x <= r
%          lb <= x <= ub
%          x belongs to BICSN

% Input objective
% ----------------
% Notice at this point the objFn is guaranteed to be scalar-valued
% if objective function is a type of variable, then input,
% otherwise, leave it as zero
prob.H = [];
add_obj = 0;
if(isa(model_obj.ObjFn, 'rome_var'))
    % zero pad to include aux model variables.
    prob.f = full(model_obj.ObjFn.BiAffineMap(1, 2:end).');
    prob.f = [zeros(model_obj.ObjFn.NumUnmappedVars, 1); prob.f];
    prob.f = [prob.f; zeros(model_obj.NumVars - length(prob.f), 1)];
    add_obj = model_obj.ObjFn.BiAffineMap(1, 1);
else
    prob.f = zeros(model_obj.NumVars, 1);
end

% check for maximimization flag
if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
    prob.f = -prob.f;
end

% Input Bounding Constraints
% --------------------------
prob.lb = [];
if(~isempty(model_obj.LB))
    prob.lb = model_obj.LB;
    num_append_LB = model_obj.NumVars - length(prob.lb);
    
    % extend constraints to fit number of variables
    if(num_append_LB > 0)
        prob.lb = [prob.lb; repmat(-Inf, num_append_LB, 1)];
    end
end

prob.ub = [];
if(~isempty(model_obj.UB))
    prob.ub = model_obj.UB;
    num_append_UB = model_obj.NumVars - length(prob.ub);
    
    % extend constraints to fit number of variables
    if(num_append_UB > 0)
        prob.ub = [prob.ub; repmat(Inf, num_append_UB, 1)];
    end
end

% Set the varriable Type (might not be needed)
% TODO: Not Implemented Yet
prob.ctype = model_obj.VarType';

% Set empty SOS constraints
prob.sos = [];

% Input Linear constraints
% ----------------
if(~isempty(model_obj.LC))
    if(~model_obj.LC.IsCertain)
        error('NOT DONE YET');
    end
    A = -model_obj.LC.BiAffineMap(:, 2:end);    % need negative sign because we have Ax + b >= 0
    A = [zeros(model_obj.LC.TotalSize, model_obj.LC.NumUnmappedVars), A];
    A = [A, zeros(size(A, 1), model_obj.NumVars - size(A, 2))]; % expand A to account for aux vars
    b = full(model_obj.LC.BiAffineMap(:, 1));   % make into full vector
else
    % If there are no linear constraints, will have to create a
    % 'fake' constraint.
    A = spalloc(1, model_obj.NumVars, 1);
    A(1, 1) = 1;
    b = Inf;
end

IndEq = false(size(A,1), 1);
IndEq(model_obj.IndEq) = true; 

prob.Aeq = A(IndEq,:);
prob.beq = b(IndEq,:);
prob.Aineq = A(~IndEq,:);
prob.bineq = b(~IndEq,:);

% Input SOC constraints
% ----------------------
% QC = [];
prob.qc = [];
if(~isempty(model_obj.QC))
    prob.qc = repmat(struct('a',[],'rhs',[],'Q',[]),1,length(model_obj.QC));

    % iterate over all the SOC constraints
    for ii = 1:length(model_obj.QC)
        cur_obj = model_obj.QC{ii};

        if(isnumeric(cur_obj))
            cur_obj = cur_obj(1):cur_obj(2);
            prob.qc(ii).Q=...
                sparse(cur_obj,cur_obj,[-1 ones(1,length(cur_obj)-1)],...
                model_obj.NumVars,model_obj.NumVars);
            prob.qc(ii).a = sparse(model_obj.NumVars,1);
            prob.qc(ii).rhs = 0;
            prob.lb(cur_obj(1)) = max(prob.lb(cur_obj(1)),0);
        else
            if(~cur_obj.IsCertain)
                error('NOT YET DONE');
            end
            
            num_pre_zeros = cur_obj.NumUnmappedVars;
            num_post_zeros= model_obj.NumVars - (cur_obj.NumMappedVars + 1 + num_pre_zeros);
            cur_A = cur_obj.BiAffineMap(2:end, 2:end);    % N-1 vars, linear component
            cur_b = cur_obj.BiAffineMap(2:end, 1);        % N-1 vars, const ccomponent
            cur_C = cur_obj.BiAffineMap(1, 2:end);        % Nth var, linear component
            cur_d = cur_obj.BiAffineMap(1, 1);            % Nth var, const component
            
            % allocate quadratic and linear terms
            prob.qc(ii).Q = spalloc(model_obj.NumVars, model_obj.NumVars, size(cur_A, 2).^2);
            prob.qc(ii).a = zeros(1, model_obj.NumVars);
            if(strcmp(cur_obj.NonLinearity, 'SQ'))
                % quad term
                prob.qc(ii).Q((num_pre_zeros+1):(end-num_post_zeros-1), ...
                    (num_pre_zeros+1):(end-num_post_zeros-1)) = cur_A' * cur_obj.ExtraData * cur_A;
                
                % linear term
                prob.qc(ii).a(num_pre_zeros + (1:cur_obj.NumMappedVars)) = 2 * (cur_b' * cur_obj.ExtraData * cur_A) - cur_C;
                
                % constant term
                prob.qc(ii).rhs = cur_d - cur_b' * cur_obj.ExtraData * cur_b ;
            else
                % ERROR
                error('rome_model:convert_cplex:UnknownNonLinearity', 'Type of Non-linearity: %s is not supported ', cur_obj.NonLinearity);
            end
        end
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
