
function [prob, add_obj]=convert_mosek(model_obj)

%ROME_MODEL\CONVERT_MOSEK Converts the model object into Mosek solvable
%form
%   Case of  LP and SOCP

%   1.Jingxi May 21


%Initialization
prob.c=[]; %objective
prob.a=[]; %linear constraints
prob.bux=[];%upper bounds on variables
prob.blx=[];%lower bounds on variables
prob.buc=[];%upper bounds on linear constraints
prob.blc=[];%lower boudns on linear constraints


if(isa(model_obj.ObjFn, 'rome_var'))
    prob.c = model_obj.ObjFn.BiAffineMap(1, 2:end).';
    prob.c = [spalloc(model_obj.ObjFn.NumUnmappedVars, 1, 0); prob.c];
    prob.c = [prob.c; spalloc(model_obj.NumVars - length(prob.c), 1,0)];  %adjust the size of prob.c to fit the number of variables 
    add_obj=model_obj.ObjFn.BiAffineMap(1,1);  %add_obj stores the constant in objective function
end

% check for maximimization flag
if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
    prob.c = -prob.c;
end

% Input Bounding Constraints
% --------------------------

%lower bound on variables
prob.blx = model_obj.LB;
% extend constraints to fit number of variables
num_append_LB = model_obj.NumVars - length(prob.blx);
if(num_append_LB > 0)
    prob.blx = [prob.blx; repmat(-Inf, num_append_LB, 1)];
end
    
%upper bound on variables
prob.bux = model_obj.UB;
% extend constraints to fit number of variables
num_append_UB = model_obj.NumVars - length(prob.bux);
if(num_append_UB > 0)
    prob.bux = [prob.bux; repmat(Inf, num_append_UB, 1)];
end
    
%linear constraints
if(~isempty(model_obj.LC))
    if(~model_obj.LC.IsCertain)
        error('NOT DONE YET');
    end
    
    %prob.a is required to be a sparse matrix
    prob.a = model_obj.LC.BiAffineMap(:, 2:end);   
    %adjust the size of prob.a to fit the number of variables
    if(model_obj.LC.NumUnmappedVars>=1)
        prob.a = [spalloc(model_obj.LC.TotalSize, model_obj.LC.NumUnmappedVars,0), prob.a];
    end
    if(model_obj.NumVars - size(prob.a, 2)>=1)
        prob.a = [prob.a, spalloc(size(prob.a, 1),model_obj.NumVars - size(prob.a, 2),0)]; 
    end
    %convert linear constraints
    prob.blc=full(-model_obj.LC.BiAffineMap(:,1)); % full added by joel 27 Mar 2010
    
    %equality constraints: b<=a'x<=b
    prob.buc=inf(size(prob.blc));
    if(~isempty(model_obj.IndEq))
%         prob.buc=inf(size(prob.blc));
        prob.buc(model_obj.IndEq,:)=prob.blc(model_obj.IndEq,:);
    end
else
    prob.a=spalloc(1,model_obj.NumVars,0);
end

 %convert quadratic constraints if any
if (~isempty(model_obj.QC))  
%     num_qc=size(model_obj.QC,1);
%     for ii=1:num_qc
%         cur_obj = model_obj.QC{ii};  % get the current SOC object
%         N = cur_obj.NumMappedVars;   % new variables to be added
%         prob.a = [prob.a, spalloc(size(prob.a, 1), N, 0)]; % extend a horz
%         D = spalloc(N, size(prob.a, 2), 2*N);
%         
%         fst_inds = cur_obj.NumUnmappedVars + (1:cur_obj.NumMappedVars);
%         D(:, fst_inds) = eye(N); 
%         
%         sec_inds = (size(prob.a, 2) - N + 1):size(prob.a, 2);
%         D(:, sec_inds) = -eye(N);
%         prob.a = [prob.a; D];
%         
%         % increase dim of c, b, etc
%         prob.blc = [prob.blc; zeros(N, 1)];
%         prob.buc = [prob.buc; zeros(N, 1)];
%         prob.blx = [prob.blx; -Inf(N, 1)];
%         prob.bux = [prob.bux;  Inf(N, 1)];
%         prob.c   = [prob.c;  zeros(N, 1)];
%                 
%         prob.cones{ii}.type='MSK_CT_QUAD';
%         prob.cones{ii}.sub=sec_inds;
%     end
    num_qc=size(model_obj.QC,1);
    for ii=1:num_qc
        prob.cones{ii}.type='MSK_CT_QUAD';
        prob.cones{ii}.sub=model_obj.QC{ii}(1):1:model_obj.QC{ii}(end);
    end
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
