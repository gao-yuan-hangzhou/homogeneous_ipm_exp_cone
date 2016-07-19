function [X_deflects, f_coeffs, bound_obj, cols] = apply_na_bdldr(model_obj, objective_obj, LDR_obj, varargin)

% ROME_MODEL\APPLY_NA_BDLDR Internal method that attempts to convert
% constraints into bi-deflected constraints, with non-anticipative
% requirements
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
% 2. This is very similar to the standard BDLDR
%
% Last Modified
% 1. Joel (9 Mar 2009)

% objective_obj defaults to model objective if none is specified
if(isempty(objective_obj))
    objective_obj = model_obj.ObjFn;
end

% check for non-linear constraints in model
if(~isempty(model_obj.QC))
    error('rome_model:apply_na_bdldr:NonLinearityPresent', ...
        'Model contains non-linearity. Deflected LDRs only work for linear models');
end

% check if a bound has been already applied
if(model_obj.bBoundUsed)
    error('rome_model:apply_na_bdldr:BoundAlreadyApplied', ... 
        'Cannot use robust bounds in conjunction with Deflected LDRs.');
end

% find indices of mean (i.e. expectation) constraints
ind_mean_constraints = model_obj.ldrMeanInd;

% extract Linear Constraints corresponding to mean constraints and standard
% LDR constraints
meanLC = model_obj.ldrLC(ind_mean_constraints, :);
ind_non_mean_constraints = setdiff(1:model_obj.ldrLC.TotalSize, ind_mean_constraints);
origLdrLC = model_obj.ldrLC(ind_non_mean_constraints, :);

% make the constraint object
constraint_obj = strip_rand(origLdrLC); % object which encodes the constraint matrix
N = length(origLdrLC);                  % total number of linear (in)equality constraints

% original index of equality constraints
origIndEq = model_obj.ldrIndEq;

% Sept 6, 2011: Joel shift index of equality constraints to account for
% mean constraints
for ii = 1:length(ind_mean_constraints)
    mod_ind = origIndEq > ind_mean_constraints(ii); % indices to be modified
    origIndEq(mod_ind) = origIndEq(mod_ind) - 1;
end

% find number of deflections (i.e. inequality constraints only)
ind_deflects = setdiff(1:N, origIndEq); % Sept 6, 2011. Joel: this is probably wrong

% obtain LDRs to deflect upon
if(isempty(LDR_obj))
    LDR_obj = model_obj.ldrVars;    
end
y0 = strip_rand(LDR_obj);

% build information index table
inner_step = LDR_obj.NumMappedVars + 1;
[r, c] = find(LDR_obj.BiAffineMap);
info_index_table = logical(accumarray(unique([r(:), ceil(c(:) ./ inner_step)], 'rows'), 1));

% construct sub-problem
y0.BiAffineMap(:, 1) = 0; % remove constant column
[rows, cols] = find(y0.BiAffineMap);
cols = unique(cols) + y0.NumUnmappedVars;   % cols counts from global 0

% sub-problem constraints
constraint_cols = cols - constraint_obj.NumUnmappedVars;
A_orig = [zeros(N, sum(constraint_cols <= 0)), ...
         -constraint_obj.BiAffineMap(:, constraint_cols(constraint_cols > 0))];        % new constraint matrix
B      = -eye(N);

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
flag_inequalities = false(N, 1);
flag_inequalities(ind_deflects) = true;
for ii = ind_deflects
    check_equal = full(sum(abs(bsxfun(@plus, A_orig(ii, :), A_orig)), 2));
    bi_deflect_table{ii} = find(check_equal == 0 & flag_inequalities);
end

% Iteration
LC_del_array = false(N, 1);             % logical array of constraints to be removed
f_coeffs = NaN(1, numel(ind_deflects)); % array of optimal values of sub problems
X_deflects = NaN(size(A_orig, 2), N); % array of optimal values of sub problems

% get verbose level and decide whether to print output
print_output = (rome_verbose > 1);

% iterate
for ii = ind_deflects    
    A = A_orig;
  
    if(~isempty(bi_deflect_table{ii}))
        % zero out unnecessary row of A
        % don't need to do the same for b, since it should be zero
        A(bi_deflect_table{ii}, :) = 0;        
    end

    % find anticipative rows and remove elements from cols
    dep_vars = (A(ii, :) ~= 0);
    cur_info = info_index_table(dep_vars(:), :); 
    cur_info = any(cur_info, 1);    % take union of dependencies
    
    % Ignore variables with NO uncertainty dependency
    if(~any(cur_info(2:end)))
        continue;
    end
    
    % find bad cols (i.e. anticipative ones) and remove them
    % from the sub-problem
    bad_cols = any(repmat(cur_info, size(info_index_table, 1), 1)...
                & (~info_index_table), 2);
            
    A(:, bad_cols) = 0;
    cur_coeff = coeff;
    cur_coeff(bad_cols) = 0;

    % check if all of A is zero
    if(all(A(:) == 0))
        continue;
    end
    
    % append equality constraint for current ii
    indEq = [origIndEq; ii];
%     [x_min, f_min, solstat, details] = cplexint([], cur_coeff, A, B(:, ii), indEq);    % solve subproblem
    switch(model_obj.Solver)
        case('CPLEXINT')     
            [x_min, f_min, solstat, details] = cplexint([], cur_coeff, A, B(:,ii), indEq);    % solve subproblem
            
            if(strcmp(details.statstring, 'optimal'))
                if(print_output)
                    disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
                end
                LC_del_array(ii) = true;
                f_coeffs(ii) = f_min;
                X_deflects(:, ii) = x_min;
            end
        % Added by Melvyn 28 May 2012:
        % For intergation CPLEX Matlab interface provided by ilog    
        case('CPLEX') 
            indInEq = setdiff((1:size(A,1))',indEq);
            [x_min, f_min, solstat,details] = cplexlp(cur_coeff, A(indInEq,:), B(indInEq,ii),A(indEq,:),B(indEq,ii));    % solve subproblem
            
            if(strcmp(details.cplexstatusstring, 'optimal'))
                if(print_output)
                    disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
                end
                LC_del_array(ii) = true;
                f_coeffs(ii) = f_min;
                X_deflects(:, ii) = x_min;
            end    
        case('MOSEK')
            prob.c=cur_coeff;
            prob.a=A;
            prob.buc=B(:,ii);
            prob.blc=ones(size(prob.buc,1),1)*-inf;
            prob.blc(indEq,1)=prob.buc(indEq,1);
            [r,res]=mosekopt('minimize echo(0)',prob);
           
            if(strcmp(res.sol.itr.solsta, 'OPTIMAL'))
                LC_del_array(ii) = true;
                f_coeffs(ii) =coeff'*res.sol.itr.xx;
                X_deflects(:, ii) =res.sol.itr.xx;
                if(print_output)
                    disp(sprintf('Constraint %d Optimal, val = %g', ii, f_coeffs(ii)));
                end                
            end
        case('SDPT3DUAL')
%             warning('Probable error using SDPT3 for deflection');
            % Debug Solve
%              [x_min, f_min, solstat, details] = cplexint([], cur_coeff, A, B(:,ii), indEq);    % solve subproblem
%             details
%              f_min
            
%             Q = eye(size(A, 1));
%             Q(:, indEq) = 0;
%             [x_min, f_min, solstat, details] = ...
%                 cplexint([], [cur_coeff; -cur_coeff; zeros(size(Q, 2), 1)], ...
%                 [A -A Q], B(:,ii), 1:size(B, 1), [], zeros(2*size(A, 2) + size(Q, 2), 1));    % solve subproblem
%             details
%             f_min
                         
%             if(strcmp(details.statstring, 'optimal'))
%                 % disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
%                 LC_del_array(ii) = true;
%                 f_coeffs(ii) = f_min;
% %                 X_deflects(:, ii) = x_min(1:size(A, 2)) + x_min(size(A, 2) + (1:size(A, 2)));
%                 xx = x_min(1:size(A, 2)) + x_min(size(A, 2) + (1:size(A, 2)));
%             end

%             % SDPT3 Primal Solve
%             Q = eye(size(A, 1));
%             Q(:, indEq) = 0;
%             At = [A, -A, Q]';
%             blk = {'l', size(At, 1)};
%             C = [cur_coeff; -cur_coeff; zeros(size(Q, 2), 1)];
%             b = B(:, ii);
% 
%             OPTIONS.printlevel=0;
%             [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);
%             if(info.termcode==0)
%                 LC_del_array(ii) = true;
%                 f_coeffs(ii) = obj(2);
%                 % reconstitute array
%                 x = X{:};
%                 X_deflects(:, ii) = x(1:length(cur_coeff)) + ... 
%                 x(length(cur_coeff) + (1:length(cur_coeff)));
%             end
% SDPT3 DUAL Solve
             % try some simple preprocessing
             [rr,cc] = find([A; cur_coeff']);
             uc = unique(cc);
             A = A(:, uc);
             cur_coeff = cur_coeff(uc);
             cur_coeff(cur_coeff == 0) = 1E-3;
             % end try
             
             NindEq = setdiff(1:N, indEq);
             blk ={'l', length(NindEq); 'u', length(indEq)};
             At  ={A(NindEq, :); A(indEq, :)};
             C   ={B(NindEq, ii); B(indEq, ii)};
             b   =-cur_coeff;
             
             OPTIONS = sqlparameters;          
             OPTIONS.printlevel=0;
             [obj,X,y,Z,info,runhist]=sdpt3(blk,At,C,b,OPTIONS);

             if(info.termcode==0)
                 LC_del_array(ii) = true;
                 f_coeffs(ii) = -obj(2);                 
                 X_deflects(:, ii) = 0;
                 X_deflects(uc, ii) = y;                
             end
             
%              % CPLEX DUAL Solve Debug
%              [x_min, f_min, solstat, details] = ...
%                  cplexint([], cat(1, C{:}), cat(1, At{:})', b, 1:length(b), [], ...
%                             [zeros(length(NindEq), 1); -inf(length(indEq), 1)]);    % solve subproblem
%              details
%              f_min
%                          
%               obj
%              keyboard
             
%             b=-cur_coeff;
%             no_constr=size(A,1);
%             rhs=B(:,ii);
%             blk{1,1}='u';
%             blk{1,2}=size(indEq,1);
%             At{1}=A(indEq,:);
%             C{1,1}=rhs(indEq,:);
%             blk{2,1}='l';
%             blk{2,2}=size(A,1)-size(indEq,1);
%             index=setdiff(1:no_constr,indEq);
%             At{2}=A(index,:);
%             C{2,1}=rhs(index,:);
            
%          OPTIONS.printlevel=0;
%         [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);

%         if(info.termcode==0)
%             LC_del_array(ii) = true;
%             f_coeffs(ii) = -obj(2);
%              X_deflects(:, ii) = y;
%             % reconstitute array
% %             x = X{:};
% %             X_deflects(:, ii) = x(1:length(cur_coeff)) +
% %             x(length(cur_coeff) + (1:length(cur_coeff)));
%         end 
    end
   
%     if(strcmp(details.statstring, 'optimal'))
%         if(print_output)
%             disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
%         end
%         LC_del_array(ii) = true;
%         f_coeffs(ii) = f_min;
%         X_deflects(:, ii) = x_min;
%     end
end

% keyboard;
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

    % TEMP
    mean_coeff(2:end, :) = -mean_coeff(2:end, :);
    opt_mean_coeff = pos(mean_coeff * X_deflects);

    % add to objective
    if(isempty(model_obj.ObjFn))
        model_obj.ObjFn = 0 ;
    end
    model_obj.ObjFn = model_obj.ObjFn + opt_mean_coeff(1, :) * bound_obj; % added ( 1, :) - big error!
    
    % add to expectation constraints.
    % Notice negative sign due to convention of NNOC versus standard
    if(~isempty(ind_mean_constraints))
        model_obj.ldrLC(ind_mean_constraints, :) = ...
            meanLC - opt_mean_coeff(1 + (1:meanLC.TotalSize), :) * bound_obj;
    end
    
    % some expectation constraints may be hidden in the certain LC.
    % Add to these guys too
    if(~isempty(ind_inequalities))
        model_obj.LC(ind_inequalities, :) = ...
            model_obj.LC(ind_inequalities, :) - opt_mean_coeff((end-numel(ind_inequalities)+1):end, :) * bound_obj;
    end

    % compute coefficients (where do I put the pos?) THINK ABOUT THIS!!
%     opt_mean_coeffs = mean_coeff * pos(X_deflects);    
% %     opt_mean_coeffs = pos(mean_coeff * X_deflects);
    
%     % add to objective
%     if(isempty(model_obj.ObjFn))
%         model_obj.ObjFn = 0 ;
%     end
%     model_obj.ObjFn = model_obj.ObjFn + opt_mean_coeffs(1, :) * bound_obj;
%     
%     % add to expectation constraints. 
%     % Notice negative sign due to convention of NNOC versus standard
%     if(~isempty(ind_mean_constraints))
%         model_obj.ldrLC(ind_mean_constraints, :) = ...
%             meanLC - neg(opt_mean_coeffs(1 + (1:meanLC.TotalSize), :)) * bound_obj;
%     end
% 
%     % some expectation constraints may be hidden in the certain LC. 
%     % Add to these guys too
%     if(~isempty(ind_inequalities))
%         model_obj.LC(ind_inequalities, :) = ...
%             model_obj.LC(ind_inequalities, :) - neg(opt_mean_coeffs((end-numel(ind_inequalities)+1):end, :)) * bound_obj;
%     end
    
    
    % shift indEq 
    newIndEqLogical = false(N, 1);
    newIndEqLogical(model_obj.ldrIndEq) = true;
    newIndEqLogical(LC_del_ind, :) = [];
    model_obj.ldrIndEq = find(newIndEqLogical);
    
    % shift ldrMeanInd
    newMeanIndLogical = false(size(model_obj.ldrMeanInd, 1) + N);
    newMeanIndLogical(model_obj.ldrMeanInd) = true;
    newMeanIndLogical(ind_non_mean_constraints(LC_del_ind), :) = [];
    model_obj.ldrMeanInd = find(newMeanIndLogical);
    
    % store results of deflection
    model_obj.ldrRemLC = model_obj.ldrLC(LC_del_ind, :);
    model_obj.ldrDefInd = cols;         % Deflection indices
    model_obj.ldrDefCoeff = X_deflects; % Deflection coefficients
    
    % remove constraints
    model_obj.ldrLC(ind_non_mean_constraints(LC_del_ind), :) = [];
    
%     %HACK
%     cols = qqq;   
%     f_coeffs = mean_coeff;
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
