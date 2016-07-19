function [x_min, f_min, sol_stat, details] = solve(model_obj, varargin)

% ROME_MODEL\SOLVE Solves the current model
%
% Assumes that at this point, all variables are included in the model
% 
% Modification History: 
% 1. Joel 
% 2. Jingxi 23 April, include MOSEK and SDPT3 as solver 
% 3. Joel 15 May 09, added default values for output_vars
% 4. Melyvn May 2012, added support for CPLEX (ilog).

 

global ROME_ENV;
x_min = 0;
f_min = 0;
sol_stat = '';
details = '';

% Override call to solve with solver choice
if(nargin > 1)
    model_obj.Solver = varargin{1};
end

% Convert into Deterministic Equivalent
convert_robust(model_obj);

% Switch based on Solver Choice
switch(model_obj.Solver)
	% MODIFIED TO CPLEXINT
    case 'CPLEXINT'
        % Primal Solve
        % format the model to conform to cplex input format
        [H, f, A, b, IndEq, QC, LB, UB, VarType, add_obj] = convert_cplexint(model_obj);
    
        % Call CPLEX to Solve
        %cplexoptions.timelimit = 120;
        %[x_min, f_min, sol_stat, details] = cplexint(H, f, A, b, IndEq, QC, LB, UB, VarType, [], cplexoptions);
        [x_min, f_min, sol_stat, details] = cplexint(H, f, A, b, IndEq, QC, LB, UB, VarType);
        
        % first check on status
        if(strcmp(details.statstring, 'unbounded or infeasible'))
            num_vars = length(VarType);
            [x_min, f_min, sol_stat, details] = cplexint(H, zeros(num_vars, 1), A, b, IndEq, QC, LB, UB, VarType);
            if(strcmp(details.statstring, 'optimal'))
                details.statstring = 'unbounded';
                x_min = NaN(num_vars, 1);
            else
                details.statstring = 'infeasible';
            end
        end
        
        % check status and flag f_min
        if(strcmp(details.statstring, 'unbounded'))
            % check for unboundedess
            f_min = -Inf;            
        elseif(strcmp(details.statstring, 'infeasible'))
            % check for infeasibility
            f_min = Inf;            
        end
        
        % Display status string
        if(ROME_ENV.VERBOSE > 0)
            disp(sprintf('Status: %s', details.statstring));
        end
        
        % if maximize flag is set, change sign of objective
        if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
            f_min = -f_min;
        end
        
        % Assign model variables
        model_obj.XSol = x_min;
        model_obj.ObjVal = f_min + add_obj;
        model_obj.DualVars = details.dual;  % notice that we don't have dualVars for x <= UB , x>=LB constraints.

	% ADDED SUPPORT FOR CPLEX MAY 2012
    case 'CPLEX'
        % Primal Solve
        % format the model to conform to cplex input format
        [prob, add_obj] = convert_cplex(model_obj);

        % Call CPLEX to Solve
        if all(prob.ctype=='C')
            [x_min,f_min,exitflag,output] = cplexqcp(prob);
        else
            [x_min,f_min,exitflag,output] = cplexmiqcp(prob);
        end

        % first check on status
        if(any(strfind(output.cplexstatusstring,'unbounded or infeasible')))
            num_vars = length(prob.f);
            prob.f = zeros(num_vars,1);
            if all(prob.ctype=='C')
                [x_min,f_min,exitflag,output] = cplexqcp(prob);
            else
                [x_min,f_min,exitflag,output] = cplexmiqcp(prob);
            end
            if(any(strfind(output.cplexstatusstring,'optimal')))
                output.cplexstatusstring = 'unbounded';
                x_min = NaN(num_vars, 1);
                f_min = -Inf; % Added by Joel June 26, 2012
            else
                output.cplexstatusstring = 'infeasible';
                f_min = Inf; % Added by Joel June 26, 2012
            end
        end
        
%         % check status and flag f_min
%         if(strcmp(output.cplexstatusstring(end-8:end),'unbounded'))
%             % check for unboundedess
%             f_min = -Inf;            
%         elseif(strcmp(output.cplexstatusstring(end-9:end),'infeasible'))
%             % check for infeasibility
%             f_min = Inf;            
%         end
        
        % Display status string
        if(ROME_ENV.VERBOSE > 0)
            disp(sprintf('Status: %s', output.cplexstatusstring));
        end
        
        % if maximize flag is set, change sign of objective
        if (model_obj.MinMaxFlag == rome_model.MAXIMIZE)
            f_min = -f_min;
        end
        
        % Assign model variables
        model_obj.XSol = x_min;
        model_obj.ObjVal = f_min + add_obj;

        % TODO: Dual vars need to be assigned.
%         model_obj.DualVars = ;  % notice that we don't have dualVars for x <= UB , x>=LB constraints.
    case 'MOSEK'
        % Convert to mosek input format
        [prob, add_obj]=convert_mosek(model_obj);

        %call MOSEK to solve
        [r,res]=mosekopt('minimize echo(0)',prob, ROME_ENV.MSK_PARAMS);
        
        % extract output
        % JOEL's NOTE(13 Sept 2010): Fails when program is unbounded. To fix this.
        if(strcmpi(res.rcodestr, 'MSK_RES_OK') || strcmpi(res.rcodestr, 'MSK_RES_TRM_STALL'))
        x_min = res.sol.itr.xx;
        f_min = prob.c'*res.sol.itr.xx;
        sol_stat = res.sol.itr.solsta;
        
        %display solution status
        if(ROME_ENV.VERBOSE > 0)
            disp(sprintf('Status: %s', res.sol.itr.solsta));
        end
        else
            % Other error in MOSEK
            disp(sprintf('MOSEK Internal Error: %s', res.rcodestr));        
        end
        
        %check for infeasibility
        if(strcmp(sol_stat, 'PRIMAL_INFEASIBLE_CER'))
            f_min = Inf;
        elseif(strcmp(sol_stat, 'DUAL_INFEASIBLE_CER'))
            f_min = -Inf;
        end
        
        if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
            f_min = -f_min;
        end
        model_obj.XSol = x_min;
        model_obj.ObjVal = f_min+add_obj;
        
    case 'SDPT3DUAL'
        % issue warning
        warning('rome_model:solve:experimental', ...
            'SDPT3DUAL Solver in Experimental Phase. Use at your own risk.');
        
        % SDPT3 dual input format
        [blk,At,C,b,add_obj]=convert_sdpt3dual(model_obj);
        
        %call SDPT3 to solve
        OPTIONS.printlevel=rome_verbose;
        [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);
        if(model_obj.MinMaxFlag == rome_model.MINIMIZE)
            obj = -obj;
        end
        model_obj.XSol=y;
        model_obj.ObjVal=obj(2)+add_obj;%add back the constant in the objective
        
        if(info.termcode==0)
            sol_stat='Optimal';
        elseif(info.termcode==1)
            sol_stat='Primal infeasible';
        elseif(info.termcode==2)
            sol_stat='Dual infeasible';
        else
            sol_stat='Unknown';
        end
        
        if(rome_verbose > 0)
            %display solution status
            disp(sprintf('\nStatus: %s', sol_stat));
        end        
        
     case 'SDPT3'
         % issue warning
         warning('rome_model:solve:experimental', ...
             'SDPT3 Solver in Experimental Phase. Use at your own risk.');
         
         %convert to SDPT3 input format
         [blk, At,C,b]=convert_sdpt3dual(model_obj);
         OPTIONS.printlevel=2;
         
         %call SDPT3 to solve
         [obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,OPTIONS);
         
         if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
             obj = -obj;
         end
         model_obj.XSol=X{1};
         model_obj.ObjVal=obj(2);
        
    case 'CPLEXDUAL'
        % issue warning
        warning('rome_model:solve:experimental', ...
            'CPLEXDUAL Solver in Experimental Phase. Use at your own risk.');
        
        % Solves primal and dual in CPLEX        
        % Primal Solve
        % format the model to conform to cplex input format
        [H, f, A, b, IndEq, QC, LB, UB, VarType] = convert_cplex(model_obj);
                      
        % Call CPLEX to Solve
        [x_min, f_min, sol_stat, details] = cplexint(H, f, A, b, IndEq, QC, LB, UB, VarType);
        
        % if maximize flag is set, change sign of objective
        if(model_obj.MinMaxFlag == rome_model.MAXIMIZE)
            f_min = -f_min;
        end
        
        % Assign model variables
        model_obj.XSol = x_min;
        model_obj.ObjVal = f_min;
        
        % Dual Solve
        % format the model to conform to cplex input format
        [H, f, A, b, IndEq, QC, LB, UB, VarType] = convert_cplexdual(model_obj);
     
        % Call CPLEX to Solve
        [x_min, f_min, sol_stat, details] = cplexint(H, f, A, b, IndEq, QC, LB, UB, VarType);
        
        % flip sign of dual objective (if primal was to minimize, then dual
        % is to maximize)
        if(model_obj.MinMaxFlag == rome_model.MINIMIZE)
            f_min = -f_min;
        end
        
        % assign model variables
        model_obj.DualVars  = x_min;
        model_obj.DualObjVal= f_min;
    otherwise
        error('rome_model:solve:InvalidSolverChoice', 'No Internal Solver Found');
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
