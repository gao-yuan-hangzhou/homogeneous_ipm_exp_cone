% function obj = rome_linearrule(sz, dep_rand_vars, varargin)
function obj = rome_linearrule(varargin)

% ROME_LINEARRULE Makes a new object of type linearrule
%
% First temp version: VERY SUBOPTIMAL - must improve this!
%
% Modified by: 
% 1. Joel

% find which argument holds the dependent uncertain variables, and assign 
dep_rand_vars = [];
numeric_flag = cellfun(@isnumeric, varargin) & ~cellfun('isempty', varargin);

% find the end of the size arguments
size_end = find(~numeric_flag, 1) - 1; 
sz = [varargin{1:size_end}];

% get the dependent variables
if(size_end < length(varargin))
    dep_rand_vars = varargin{size_end+1};
end

% get the options
varargin = varargin(size_end+2:end);

% check if options contain the 'Pattern' string
LDRpattern = [];

% do this check if varargin is non-empty
if(~isempty(varargin))
    patternFlag = cellfun(@strcmp, varargin, cellstr(repmat('Pattern', numel(varargin), 1))', 'UniformOutput', true);
    if(any(patternFlag))
        % check that only 1 pattern can be specified
        if(sum(patternFlag) > 1)
            error('rome_linearrule:OnlyOnePattern', 'Cannot specify more than 1 pattern');
        end
        % check that pattern is logical
        if(~islogical(patternFlag))
            error('rome_linearrule:LogicalArrayReq', 'Pattern has to be a logical array');
        end
        patternInd = find(patternFlag);
        LDRpattern = varargin{patternInd+1};
        varargin(patternInd + [0, 1]) = [];
    end
end
% default size values
if(isempty(sz))
    sz = [1, 1];
elseif(numel(sz) < 2)
    sz = [sz, 1];
end

if(isempty(dep_rand_vars))
    obj = rome_model_var(sz);
else
    if(~dep_rand_vars.IsRand)
        error('rome_linearrule:MustBeUncertain', 'Dependent Uncertain Vars must be pure uncertain!');
    end

    % vectorize randvars
    dep_rand_vars = dep_rand_vars(:);   % vectorize

    % optimization for diagonal rand_vars
    if(dep_rand_vars.IsDiag)
        % number of cols in each block
        inner_sz = prod(sz) * (dep_rand_vars.TotalSize + 1);

        obj = rome_model_var(inner_sz);
        obj.Size = sz;
        obj.NumUnmappedRandVars = dep_rand_vars.NumUnmappedRandVars;
        obj.NumMappedRandVars = dep_rand_vars.TotalSize;
        
        % manually define BiAffineMap
        obj.BiAffineMap = [zeros(inner_sz, 1), obj.DiagMult .* speye(inner_sz)];
        obj.BiAffineMap = blk_straighten(obj.BiAffineMap, prod(sz), inner_sz + 1, 'row');
        
        % apply pattern if non-empty
        if(~isempty(LDRpattern))
            % check that sizing is correct
            if((size(LDRpattern, 1) ~= prod(sz)) || (size(LDRpattern, 2) ~= (1 + length(dep_rand_vars))))
                error('rome_linearrule:InvalidPatternSize', 'Pattern must be a 2D matrix with dimensions TotalSize x (1 + NDependentVars)');
            end

            % Make an expanded pattern
            LDRpatternExpand = repmat(LDRpattern(:), 1, obj.NumMappedVars + 1);
            LDRpatternExpand = blk_transpose(LDRpatternExpand, size(LDRpattern, 1), obj.NumMappedVars + 1);

            % apply pattern to the BiAffineMap
            obj.BiAffineMap(~LDRpatternExpand) = 0;
        end
    else

        % apply pattern if non-empty
        if(~isempty(LDRpattern))
            % check that sizing is correct
            if((size(LDRpattern, 1) ~= prod(sz)) || (size(LDRpattern, 2) ~= (1 + length(dep_rand_vars))))
                error('rome_linearrule:InvalidPatternSize', 'Pattern must be a 2D matrix with dimensions TotalSize x (1 + NDependentVars)');
            end

            % allocate space for object
            obj = rome_model_var(sz);
            obj.BiAffineMap(~LDRpattern(:, 1), :) = 0;
            
            % standard rand vars (no pattern)
            for ii = 1:length(dep_rand_vars)
                new_obj = rome_model_var(sz);
                new_obj.BiAffineMap(~LDRpattern(:, 1+ii), :) = 0;
                obj = obj + dep_rand_vars(ii) .* new_obj;
            end
        else
            % allocate space for object
            obj = rome_model_var(sz);
            
            % standard rand vars (no pattern)
            for ii = 1:length(dep_rand_vars)
                obj = obj + dep_rand_vars(ii) .* rome_model_var(sz);
            end
        end
    end
    
    % register object with model
    h = rome_get_current_model();
    h.ldrVars = [h.ldrVars; obj(:)]; 
end

% parse optional arguments
for ii = 1:2:numel(varargin)
    if(~ischar(varargin{ii}) || ~isnumeric(varargin{ii+1}))
        error('rome_linearrule:InvalidOps', 'Optional Arguments must be entered in pairs');
    else
        eval(sprintf('obj.%s = %g;', varargin{ii}, varargin{ii+1}));
    end
end


% apply constraint
rome_constraint(obj);




% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
