% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.function apply_bdldr(model_obj, objective_obj, LDR_obj, varargin)

% ROME_MODEL\APPLY_BDLDR Internal method that attempts to convert
% constraints into bi-deflected constraints
%
% Inputs:
%   model_obj      : Current model object
%   objective_obj  : A rome_var object representing the objective function
%   LDR_obj        : An LDR rome_var object represting the variables to be
%                    deflected
%   varargin       : A list of function handles of bounds to be used to
%                    approximate E()^+
%
% Notes:
% 1. The reason for the rather odd construction of the sub-problem
%    constraints and objective, is that we need to zero-out the unaccounted
%    variables (i.e. variables that are unmapped by the constraint object
%    or the objective object).
%
% 2. This is very similar to the DLDR case.
%
% Last Modified
% 1. Joel (9 Mar 2009)

% find indices of mean (i.e. expectation) constraints
ind_mean_constraints = model_obj.ldrMeanInd;

% extract Linear Constraints corresponding to mean constraints and standard
% LDR constraints
meanLC = model_obj.ldrLC(ind_mean_constraints, :);
origLdrLC = model_obj.ldrLC(setdiff(1:model_obj.ldrLC.TotalSize, ind_mean_constraints), :);

% make the constraint object
constraint_obj = strip_rand(origLdrLC); % object which encodes the constraint matrix
N = length(origLdrLC);                  % total number of linear (in)equality constraints

% original index of equality constraints
origIndEq = model_obj.ldrIndEq;

% find number of deflections (i.e. inequality constraints only)
ind_deflects = setdiff(1:N, origIndEq);

% obtain LDRs to deflect upon
if(isempty(LDR_obj))
    y0 = strip_rand(origLdrLC);    % LDRs to deflect on
else
    y0 = strip_rand(LDR_obj);      % LDRs to deflect on
end

% construct sub-problem
y0.BiAffineMap(:, 1) = 0; % remove constant column
[rows, cols] = find(y0.BiAffineMap);
cols = unique(cols) + y0.NumUnmappedVars;   % cols counts from global 0

% sub-problem constraints
constraint_cols = cols - constraint_obj.NumUnmappedVars;
A_orig = [zeros(N, sum(constraint_cols <= 0)), ...
         -constraint_obj.BiAffineMap(:, constraint_cols(constraint_cols > 0))];        % new constraint matrix
B      = -eye(N);

% construct sub-problem
y0.BiAffineMap(:, 1) = 0; % remove constant column
[rows, cols] = find(y0.BiAffineMap);
cols = unique(cols) + y0.NumUnmappedVars;   % cols counts from global 0

% sub-problem objective
if(isnumeric(objective_obj) || isempty(objective_obj))
    coeff = zeros(size(A_orig, 2), 1);
else
    certain_objective = strip_rand(objective_obj);
    certain_objective.BiAffineMap(:, 1) = 0; % remove constant column

    objective_cols = cols - certain_objective.NumUnmappedVars;
    valid_cols = objective_cols(objective_cols > 0 & objective_cols <= certain_objective.NumMappedVars + 1);
    coeff = [zeros(1, sum(objective_cols <= 0)),  ...
        certain_objective.BiAffineMap(:, valid_cols), ...
        zeros(1, sum(objective_cols > certain_objective.NumMappedVars + 1))]; % coefficients of objective function
    coeff = full(coeff'); % vectorize
end

% handle mean-constraints
ind_inequalities = [];
if(~isempty(model_obj.LC))
    ind_inequalities = setdiff(1:model_obj.LC.TotalSize, model_obj.IndEq);
end
certain_mean = strip_rand(meanLC);
certain_mean = [model_obj.ObjFn; certain_mean; model_obj.LC(ind_inequalities, :)]; % prepend model objective
certain_mean.BiAffineMap(:, 1) = 0; % remove constant column


mean_cols = cols - certain_mean.NumUnmappedVars;
valid_mean_cols = mean_cols(mean_cols > 0 & mean_cols <= certain_mean.NumMappedVars + 1);
mean_coeff = [zeros(certain_mean.TotalSize, sum(mean_cols <= 0)),  ...
         certain_mean.BiAffineMap(:, valid_mean_cols), ...
         zeros(certain_mean.TotalSize, sum(mean_cols > certain_mean.NumMappedVars + 1))]; % coefficients of mean constraints
mean_coeff = full(mean_coeff); % notice no transpose

% find bi-deflections
bi_deflect_table = cell(N, 1);
for ii = ind_deflects
    check_equal = sum(abs(bsxfun(@plus, A_orig(ii, :), A_orig)), 2);
    bi_deflect_table{ii} = find(check_equal == 0);
end

% Iteration
LC_del_array = false(N, 1);             % logical array of constraints to be removed
f_coeffs = NaN(1, numel(ind_deflects)); % array of optimal values of sub problems
X_deflects = NaN(size(A_orig, 2), numel(ind_deflects)); % array of optimal values of sub problems
for ii = ind_deflects
    A = A_orig;
    
    if(~isempty(bi_deflect_table{ii}))
        % zero out unnecessary row of A
        % don't need to do the same for b, since it should be zero
        A(bi_deflect_table{ii}, :) = 0;        
    end    
    
    indEq = [origIndEq; ii];
     switch(model_obj.Solver)
         case('CPLEXINT')
             [x_min, f_min, solstat, details] = cplexint([], coeff, A, B(:, ii), indEq);    % solve subproblem
             
             if(strcmp(details.statstring, 'optimal'))
                 disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
                 LC_del_array(ii) = true;
                 f_coeffs(ii) = f_min;
                 X_deflects(:, ii) = x_min;
             end
        % Added by Melvyn 28 May 2012:
        % For intergation CPLEX Matlab interface provided by ilog    
        case('CPLEX') 
            indInEq = setdiff((1:size(A,1))',indEq);
            [x_min, f_min, solstat,details] = cplexlp(coeff, A(indInEq,:), B(indInEq,ii),A(indEq,:),B(indEq,ii));    % solve subproblem
            
            if(strcmp(details.cplexstatusstring, 'optimal'))
                if(print_output)
                    disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
                end
                LC_del_array(ii) = true;
                f_coeffs(ii) = f_min;
                X_deflects(:, ii) = x_min;
            end    
          case('MOSEK')
            prob.c=coeff;
            prob.a=A;
            prob.buc=B(:,ii);
            prob.blc=ones(size(prob.buc,1),1)*-inf;
            prob.blc(indEq,1)=prob.buc(indEq,1);
            [r,res]=mosekopt('minimize echo(0)',prob);
           if(strcmp(res.sol.itr.solsta, 'OPTIMAL'))
            LC_del_array(ii) = true;
            f_coeffs(ii) =coeff'*res.sol.itr.xx;
            X_deflects(:, ii) =res.sol.itr.xx;
           disp(sprintf('Constraint %d Optimal, val = %g', ii,f_coeffs(ii)));
           end 
         case('SDPT3DUAL')
            b=-coeff;
            no_constr=size(A,1);
            rhs=B(:,ii);
            blk{1,1}='u';
            blk{1,2}=size(indEq,1);
            At{1}=A(indEq,:);
            C{1,1}=rhs(indEq,:);
            blk{2,1}='l';
            blk{2,2}=size(A,1)-size(indEq,1);
            index=setdiff(1:no_constr,indEq);
            At{2}=A(index,:);
            C{2,1}=rhs(index,:);
           
         OPTIONS.printlevel=0;
        [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);
        if(info.termcode==0)
            LC_del_array(ii) = true;
            f_coeffs(ii) =-obj(2);
            X_deflects(:, ii) =y;
        end 
             
     end
end

% remove infeasible points
f_coeffs(isnan(f_coeffs)) = [];
X_deflects(:, isnan(X_deflects(1, :)))  = [];

if(isempty(f_coeffs))    
    disp('No Feasible Deflections.');
else
    % only run here if there was at least one feasible 
    % convert to indices
    LC_del_ind = find(LC_del_array);

    % make bound
    r = origLdrLC(LC_del_ind, :);    
    bound_obj = rome_create_bound(-r, varargin{:});

    % add to expectation constraints. 
    % Notice negative sign due to convention of NNOC versus standard
    opt_mean_coeffs = mean_coeff * pos(X_deflects); 
    
    % add to objective
    model_obj.ObjFn = model_obj.ObjFn + opt_mean_coeffs(1, :) * bound_obj ;
    
    % add to expectation constraints. 
    % Notice negative sign due to convention of NNOC versus standard
    if(~isempty(ind_mean_constraints))
        model_obj.ldrLC(ind_mean_constraints, :) = ...
            meanLC - neg(opt_mean_coeffs(1 + (1:meanLC.TotalSize), :)) * bound_obj;
    end

    % some expectation constraints may be hidden in the certain LC. 
    % Add to these guys too
    if(~isempty(ind_inequalities))
        model_obj.LC(ind_inequalities, :) = ...
            model_obj.LC(ind_inequalities, :) - neg(opt_mean_coeffs((end-numel(ind_inequalities)+1):end, :)) * bound_obj;
    end    
    
    % shift indEq 
    newIndEqLogical = false(N, 1);
    newIndEqLogical(model_obj.ldrIndEq) = true;
    newIndEqLogical(LC_del_ind, :) = [];
    model_obj.ldrIndEq = find(newIndEqLogical);
    
    % store results of deflection
    model_obj.ldrRemLC = model_obj.ldrLC(LC_del_ind, :);
    model_obj.ldrDefInd = cols;         % Deflection indices
    model_obj.ldrDefCoeff = X_deflects; % Deflection coefficients
    
    % remove constraints
    model_obj.ldrLC(LC_del_ind, :) = [];
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
