function [X_deflects, f_coeffs, cols] = apply_dldr(model_obj, objective_obj, LDR_obj, varargin)

% ROME_MODEL\APPLY_DLDR Internal method that attempts to convert
% constraints into deflected constraints
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
% Last Modified
% 1. Joel (9 Mar 2009)

% get objects
origLdrLC = model_obj.ldrLC;            % object of Linear constraints
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

% get a certain version of the objective object
% certain_objective_obj = strip_rand(objective_obj);

% construct sub-problem
y0.BiAffineMap(:, 1) = 0; % remove constant column
[rows, cols] = find(y0.BiAffineMap);
cols = unique(cols) + y0.NumUnmappedVars;   % cols counts from global 0

% sub-problem constraints
constraint_cols = cols - constraint_obj.NumUnmappedVars;
A_orig = [zeros(constraint_obj.TotalSize, sum(constraint_cols <= 0)), ...
         -constraint_obj.BiAffineMap(:, constraint_cols(constraint_cols > 0))];        % new constraint matrix
B      = -eye(N);

% sub-problem objective
certain_objective = strip_rand(objective_obj);
certain_objective.BiAffineMap(:, 1) = 0; % remove constant column

objective_cols = cols - certain_objective.NumUnmappedVars;
valid_cols = objective_cols(objective_cols > 0 & objective_cols <= certain_objective.NumMappedVars + 1);
coeff = [zeros(1, sum(objective_cols <= 0)),  ...
         certain_objective.BiAffineMap(:, valid_cols), ...
         zeros(1, sum(objective_cols > certain_objective.NumMappedVars + 1))]; % coefficients of objective function
coeff = full(coeff'); % vectorize

% Iteration
LC_del_array = false(N, 1);             % logical array of constraints to be removed
f_coeffs = NaN(1, numel(ind_deflects)); % array of optimal values of sub problems
X_deflects = NaN(size(A_orig, 2), numel(ind_deflects)); % array of optimal values of sub problems
for ii = ind_deflects
    indEq = [origIndEq; ii];
    switch(model_obj.Solver)
        case('CPLEXINT')
    [x_min, f_min, solstat, details] = cplexint([], coeff, A_orig, B(:,ii), indEq);    % solve subproblem
    
    if(strcmp(details.statstring, 'optimal'))
        % disp(sprintf('Constraint %d Optimal, val = %g', ii, f_min));
        LC_del_array(ii) = true;
        f_coeffs(ii) = f_min;
        X_deflects(:, ii) = x_min;
    end
        % Added by Melvyn 28 May 2012:
        % For intergation CPLEX Matlab interface provided by ilog    
        case('CPLEX') 
            indInEq = setdiff((1:size(A_orig,1))',indEq);
            [x_min, f_min, solstat,details] = cplexlp(coeff, A_orig(indInEq,:), B(indInEq,ii),A_orig(indEq,:),B(indEq,ii));    % solve subproblem
            
            if(strcmp(details.cplexstatusstring, 'optimal'))
                LC_del_array(ii) = true;
                f_coeffs(ii) = f_min;
                X_deflects(:, ii) = x_min;
            end    
        case('MOSEK')
            prob.c=coeff;
            prob.a=A_orig;
            prob.buc=B(:,ii);
            prob.blc=ones(size(prob.buc,1),1)*-inf;
            prob.blc(indEq,1)=prob.buc(indEq,1);
            [r,res]=mosekopt('minimize echo(0)',prob);
           if(strcmp(res.sol.itr.solsta, 'OPTIMAL'))
            LC_del_array(ii) = true;
            f_coeffs(ii) =coeff'*res.sol.itr.xx;
            X_deflects(:, ii) =res.sol.itr.xx;
           end 
        case('SDPT3DUAL')
            A = A_orig;
            cur_coeff = coeff;

            % try some simple preprocessing            
            [rr,cc] = find([A; cur_coeff']);
            uc = unique(cc);
            A = A(:, uc);
            cur_coeff = cur_coeff(uc);            
            cur_coeff(cur_coeff == 0) = 1E-3;
            
            [rr,cc] = find([A, B(:, ii)]);
            ur = unique(rr);
            A = A(ur, :);
            rhs = B(ur, ii);
            
            [rr, cc] = find(A);
            zrA = setdiff(size(A, 1), rr); % find rows of A with all zeros
            if(any(rhs(zrA) < 0))
                % immediately infeasible
                continue;
            end            
            indEq = intersect(ur, indEq);
            NindEq = setdiff(ur, indEq);            
            % end try
            
            blk ={'l', length(NindEq); 'u', length(indEq)};
            At  ={A(NindEq, :); A(indEq, :)};
            C   ={rhs(NindEq); rhs(indEq)};
            b   =-cur_coeff;
            
            OPTIONS = sqlparameters;
            OPTIONS.printlevel=0;
            [obj,X,y,Z,info,runhist]=sdpt3(blk,At,C,b,OPTIONS);
            
            if(info.termcode==0)
%                 keyboard
                LC_del_array(ii) = true;
                f_coeffs(ii) = -obj(2);
                X_deflects(:, ii) = 0;
                X_deflects(uc, ii) = y;
            end
%             b=-coeff;
%             no_constr=size(A_orig,1);
%             rhs=B(:,ii);
%             blk{1,1}='u';
%             blk{1,2}=size(indEq,1);
%             At{1}=A_orig(indEq,:);
%             C{1,1}=rhs(indEq,:);
%             blk{2,1}='l';
%             blk{2,2}=size(A_orig,1)-size(indEq,1);
%             index=setdiff(1:no_constr,indEq);
%             At{2}=A_orig(index,:);
%             C{2,1}=rhs(index,:);
%            
%          OPTIONS.printlevel=0;
%         [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);
%         if(info.termcode==0)
%             LC_del_array(ii) = true;
%             f_coeffs(ii) =-obj(2);
%             X_deflects(:, ii) =y;
%         end 
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
    add_obj = f_coeffs * rome_create_bound(-r, varargin{:});

    % add to objective
    model_obj.ObjFn = model_obj.ObjFn + add_obj;

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
