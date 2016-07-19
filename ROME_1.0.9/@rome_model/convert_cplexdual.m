function [H, f, A, b, IndEq, QC, LB, UB, VarType] = convert_cplexdual(model_obj)

% ROME_MODEL\CONVERT_CPLEXDUAL Converts format from model object into dual
% solvable form
% 
% Modification History: 
% 1. Joel 


% format the model in dual form to conform to cplex input format
H = [];

% find non-inf bounds
primal_LB_ind = find(~isinf(model_obj.LB));
primal_UB_ind = find(~isinf(model_obj.UB));
num_LB = length(primal_LB_ind);
num_UB = length(primal_UB_ind);

% find number of primal constraints
num_primal_LC = size(model_obj.LC.AffineMap, 1);

% number of dual variables
num_dual_vars = num_primal_LC + num_LB + num_UB;

% Input objective
% ----------------
% trying to maximize dual, so double negation
f = -full([-model_obj.LC.AffineMap(:, end); model_obj.LB(primal_LB_ind); model_obj.UB(primal_UB_ind)]);

% Input Linear constraints
% ----------------
if(~isempty(model_obj.LC))
    % initialize constraint matrix
    A = model_obj.LC.AffineMap(:, 1:end-1).';

    % number of dual linear constraints
    num_dual_LC = size(A, 1);

    % make selection matrices for upper and lower bounds
    if(num_LB ~= 0)
        LB_sel_mat = spalloc(num_dual_LC, num_LB, num_LB);
        LB_sel_mat(sub2ind(size(LB_sel_mat), primal_LB_ind.', 1:num_LB)) = 1;
    else
        LB_sel_mat = [];
    end
    if(num_UB ~= 0)
        UB_sel_mat = spalloc(num_dual_LC, num_UB, num_UB);
        UB_sel_mat(sub2ind(size(UB_sel_mat), primal_UB_ind.', 1:num_UB)) = 1;
    else
        UB_sel_mat = [];
    end

    % put into final linear constraint matrix
    A = [A, LB_sel_mat, UB_sel_mat];

    % make linear constraint coefficient
    b = full(model_obj.ObjFn.AffineMap(1, 1:end-1)).';

    % if primal was a maximization problem, this must be flipped
    if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
        b = -b;
    end
end

% Define equalities
% Joel: probably can do some optimization here to remove unncessary
% constraints
IndEq = (1:model_obj.NumVars).';

% Input Bounding Constraints
% --------------------------
LB = zeros(num_dual_vars, 1);
LB(model_obj.IndEq) = -Inf;            % Unbounded vars for Primal Equalities
LB((end - num_UB + 1): end) = -Inf;    % No LBs for lambda_U

UB = Inf(num_dual_vars, 1);            % most vars have no UB
UB((end - num_UB + 1): end) = 0;       % set UB for lambda_U

% Set the variable Type (might not be needed)
% only continous variables allowed here for dual
VarType = repmat('C', num_dual_vars, 1);

% Input SOC constraints
% ----------------------
QC = [];


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
