function convert_robust(model_obj)

% ROME_MODEL\CONVERT_ROBUST Converts the robust parameters in the model into
% standard deterministic parameters
% 
% Modification History: 
% 1. Joel 
% 2. Jingxi 14 April 09, line 45, line 196


% dbstop in convert_robust.m at 67

% check if empty
if(isempty(model_obj.ldrLC)) 
    return;
end

% expand rndUB / LB arrays
rndLB = model_obj.rndLB;
num_append_LB = model_obj.NumRandVars - length(rndLB);
rndLB = [rndLB; -Inf(num_append_LB * (num_append_LB > 0), 1)]; % extend constraints to fit number of variables

rndUB = model_obj.rndUB;
num_append_UB = model_obj.NumRandVars - length(rndUB);
rndUB = [rndUB;  Inf(num_append_UB * (num_append_UB > 0), 1)]; % extend constraints to fit number of variables

% for time being, don't allow infinite constraints
% if(any(isinf(rndLB)) || any(isinf(rndUB)))
%     error('NOT YET DONE');
% end

% Deal with Equality Constraints
% ------------------------------
if(~isempty(model_obj.ldrIndEq))
    % overload subscripted ref operator and select out LDR Equalities
    ldrEqLC = model_obj.ldrLC(model_obj.ldrIndEq, :); 
    
    % define for convenience
    block_width = ldrEqLC.NumMappedVars + 1;
    block_height = ldrEqLC.TotalSize;
    
    % Joel: Modified to handle general equalities on the uncertainties
    if(~isempty(model_obj.rndIndEq))
        % get the uncertain equality constraints
        rndEq = model_obj.rndLC(model_obj.rndIndEq);
        
        % find null space of the linear part of the affinemap
        full_linearmap = [zeros(rndEq.TotalSize, pos(ldrEqLC.NumMappedRandVars - rndEq.NumMappedRandVars)), ...
                            full(rndEq.BiAffineMap(:, 2:end))];
        nullRndEq = null(full_linearmap);
        
        % find a particular solution
        z_particular = -full_linearmap \ rndEq.BiAffineMap(:, 1);
        
        % combine and truncate
        z_new = [z_particular, nullRndEq];
        start_ind = ldrEqLC.NumUnmappedRandVars + 1;
        end_ind   = start_ind + ldrEqLC.NumMappedRandVars - 1;
        z_new = [ones(1, size(z_new, 2)); z_new(start_ind:end_ind, :)];
        
        % reshape LdrLCEq coefficient matrix by stretching each block into
        % a column
        reshaped_coeff_matrix = reshape(ldrEqLC.BiAffineMap, block_height * block_width, 1 + ldrEqLC.NumMappedRandVars);
        
        % create new coefficient matrix by multiplying (it's really a
        % projection into a lower dimensional space)
        new_coeff_matrix = reshaped_coeff_matrix * z_new;
        
        % now unshape and set
        ldrEqLC.BiAffineMap = reshape(new_coeff_matrix, block_height, block_width * size(z_new, 2));
    end

    
%     if (~isempty(model_obj.rndEq))
%         % Remove any 'false' uncertain variables
% 
%         %Jingxi :first check if 'false' uncertain variables are mapped in ldrEqLC
%         begin_ind=find(model_obj.rndEq(:,1)>ldrEqLC.NumUnmappedRandVars,1);   % Jingxi: starts from the first mapped 'false' rand variable
%         if (~isempty(begin_ind))
%             end_ind=find(model_obj.rndEq(:,1)>ldrEqLC.NumUnmappedRandVars+size(ldrEqLC, 1),1)-1;
% 
%             if (isempty(end_ind))
%                 end_ind=size(model_obj.rndEq, 1);   %Jingxi:ends at the last mapped 'false' rand variable
%             end
% 
%             %for ii = 1:size(model_obj.rndEq, 1)
%             %
%             for ii=begin_ind:end_ind
%                 % get index of block
%                 cur_block_ind = model_obj.rndEq(ii, 1) + 1 - ldrEqLC.NumUnmappedRandVars;
%                 % multiply and accumulate into first block
%                 ldrEqLC.BiAffineMap(:, 1:block_width) = ldrEqLC.BiAffineMap(:, 1:block_width) + ...
%                     model_obj.rndEq(ii, 2) * ldrEqLC.BiAffineMap(:, (block_width * (cur_block_ind-1) + 1) : block_width * cur_block_ind);
% 
%                 % zero out false block
%                 ldrEqLC.BiAffineMap(:, (block_width * (cur_block_ind-1) + 1) : block_width * cur_block_ind) = ...
%                     spalloc(block_height, block_width, 0);
% 
%             end
%         end
%     end
    
    % reshape the BiAffine Map to make new variables
    ldrEqLC.BiAffineMap = blk_transpose(ldrEqLC.BiAffineMap, block_height, block_width);
    
    % remove all rows with only zeros
    row_sum_abs = sum(abs(ldrEqLC.BiAffineMap), 2);
    ldrEqLC.BiAffineMap = ldrEqLC.BiAffineMap(row_sum_abs ~= 0, :);
    
    % Manually change the size of the rome_var to reflect the change
    ldrEqLC.Size = size(ldrEqLC.BiAffineMap, 1);
%     ldrEqLC.Size = [ldrEqLC.TotalSize * (ldrEqLC.NumMappedRandVars + 1), 1];
    
 
    % Remove all Mapped Uncertain Vars
    ldrEqLC.NumMappedRandVars = 0;
    
    % Set the cone to zero
    ldrEqLC.Cone = rome_constants.ZERO;
    
    % apply constraint 
    rome_constraint(ldrEqLC);
    
    
end

% remove equality constraints
model_obj.ldrLC = model_obj.ldrLC(setdiff(1:size(model_obj.ldrLC, 1), model_obj.ldrIndEq));
model_obj.ldrIndEq = [];

% Convert Inequality Constraints
% ------------------------------
% check if empty
if(isempty(model_obj.ldrLC.BiAffineMap)) 
    return;
end

% object representing the current constraint in the loop
obj = model_obj.ldrLC; %(ii, :);

% define params for convenience
N = obj.NumMappedRandVars;

% At this point, might want to remove rows of zeros to eliminate
bounds_ind = 1:model_obj.NumRandVars;

% get actual lower bounds
actualLB = rndLB(bounds_ind);
actualUB = rndUB(bounds_ind);

% finite bounds
indLBfinite = find(~isinf(actualLB));
indUBfinite = find(~isinf(actualUB));

% define default values
LCterm    = 0;
LCEqterm  = 0;
LCAterm   = 0;
LCAeqterm = 0;

% find (pure) uncertain linear constraints
if(~isempty(model_obj.rndLC))
    rndIndEq = false(model_obj.rndLC.TotalSize, 1);
    rndIndEq(model_obj.rndIndEq) = true;

    % get the multipliers
    rndAeq = model_obj.rndLC.BiAffineMap(rndIndEq, 2:end);
    rndbeq = -model_obj.rndLC.BiAffineMap(rndIndEq, 1); 
    rndA = model_obj.rndLC.BiAffineMap(~rndIndEq, 2:end);
    rndb = -model_obj.rndLC.BiAffineMap(~rndIndEq, 1);     
    
    % Define dual terms
    if(~isempty(rndb))
        % attach prezeros and postzeros
        num_pre_zeros = model_obj.rndLC.NumUnmappedRandVars;
        num_post_zeros = model_obj.NumRandVars - (num_pre_zeros + model_obj.rndLC.NumMappedRandVars);
        rndA = [spalloc(size(rndA, 1), num_pre_zeros, 0), ...
                rndA, ...
                spalloc(size(rndA, 1), num_post_zeros, 0)];
        
        lambda = rome_model_var(length(rndb), obj.TotalSize, 'Cone', rome_constants.NNOC);
        LCterm    = lambda' * rndb;
        LCAterm   = rndA'   * lambda;
    end
    
    if(~isempty(rndbeq))
        % attach prezeros and postzeros
        num_pre_zeros = model_obj.rndLC.NumUnmappedRandVars;
        num_post_zeros = model_obj.NumRandVars - (num_pre_zeros + model_obj.rndLC.NumMappedRandVars);
        rndAeq = [spalloc(size(rndAeq, 1), num_pre_zeros, 0), ...
                   rndAeq, ...
                   spalloc(size(rndAeq, 1), num_post_zeros, 0)];
            
        mu     = rome_model_var(length(rndbeq), obj.TotalSize);
        LCEqterm  = mu'     * rndbeq;
        LCAeqterm = rndAeq' * mu;
    end 
end

QCterm = 0;
QCAterm = 0;
if(~isempty(model_obj.rndQC))
    for ii = 1:length(model_obj.rndQC)
        cur_obj = model_obj.rndQC{ii};
        rndA = cur_obj.BiAffineMap(:, 2:end);
        rndb = -cur_obj.BiAffineMap(:, 1);
%         lambda = rome_model_var(length(rndb), 'Cone', rome_constants.SOC);
        lambda = rome_model_var(length(rndb), obj.TotalSize, 'Cone', rome_constants.SOC);
        
        num_pre_zeros = cur_obj.NumUnmappedRandVars;
        num_post_zeros = model_obj.NumRandVars - (num_pre_zeros + cur_obj.NumMappedRandVars);
        
        % attach prezeros and postzeros
        rndA = [spalloc(size(rndA, 1), num_pre_zeros, 0), ...
                rndA, ...
                spalloc(size(rndA, 1), num_post_zeros, 0)];
        
        QCAterm = QCAterm + rndA' * lambda;
        QCterm = QCterm + lambda' * rndb;
    end    
end

% define multiplier variables for lb and ub
num_dual_vars = model_obj.NumRandVars;
r = rome_empty_var(num_dual_vars, obj.TotalSize);
r(indLBfinite, :) = rome_model_var(numel(indLBfinite), obj.TotalSize, 'Cone', rome_constants.NNOC);

s = rome_empty_var(num_dual_vars, obj.TotalSize);
s(indUBfinite, :) = rome_model_var(numel(indUBfinite), obj.TotalSize, 'Cone', rome_constants.NNOC);

% Extract x0 and x_vec components from the constraint object
x_0   = strip_rand(obj);
x_vec = strip_certain(obj);

% perform inner product on finite bounds only
LBterm   = r(indLBfinite, :)' * actualLB(indLBfinite);
UBterm   = s(indUBfinite, :)' * actualUB(indUBfinite);
rome_constraint(LBterm - UBterm + LCterm + LCEqterm + QCterm + x_0 >= 0);

% add a new constraint
%Zblank = spalloc(num_dual_vars - N,obj.TotalSize, 0);

Zblank_pre=spalloc(obj.NumUnmappedRandVars,obj.TotalSize,0);% Jingxi: zero-pad for unmapped rand variables before obj
Zblank_post=spalloc(pos(num_dual_vars-N-obj.NumUnmappedRandVars),obj.TotalSize,0);% Jingxi: zero-pad for unmapped rand variables after obj

%rome_constraint((r - s) + LCAterm + LCAeqterm == [x_vec; Zblank]);
% test
% if(~isnumeric(QCAterm))
%     QCAterm = repmat(QCAterm, [1, obj.TotalSize]);
% end
rome_constraint((r - s) + LCAterm + LCAeqterm + QCAterm == [Zblank_pre;x_vec; Zblank_post]);

% remove all ldr-type constraints
model_obj.ldrLC = [];

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
