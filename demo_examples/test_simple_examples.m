addpath(fileparts(pwd)); addpath([fileparts(pwd), '/subroutines']);
while(1)
% Example 1
% Covnert a log optimization problem into a standard conic program
% and solve it using hsd_lqe
% consider the following problem (P)
%       min (-log(x1) - 2*log(x2)) 
%       s.t. x1+x2=5, x1>=0, x2>=0
% which has an optimal solution x1=5/3, x2=10/3 
% and optimal objective -2.9188 by basic calculus.

% The problem (P) can be converted into a standard conic program (P')
%       min -t(1) - 2t(2)
%       s.t. A[x_bar_1, x_bar_2]=b
% where A = [0,1,0,0,1,0; 
%            0,0,1,0,0,0; 
%            0,0,0,0,0,1]
% x_bar_1 = [t1;x1;x3] in K_exp, x_bar_2 = [t2;x2;x4], b = [5;1;1] in K_exp

% Set the input arguments for hsd_lqe
clear blk At c;
blk{1,1} = 'e'; blk{1,2} = 3*ones(2,1); At{1} = sparse([0,1,0,0,1,0; 0,0,1,0,0,0; 0,0,0,0,0,1]); c{1} = [-1;0;0;-2;0;0]; 
b = [5;1;1];
% Solve the problem using hsd_lqe
[opt_obj, x_return, y_return, z_return, info] = hsd_lqeu_fast(blk, At, c,b);
display(['x1_opt = ' num2str(x_return{1}(2)), ', x2_opt = ' num2str(x_return{1}(5))]);

m_choice=input('Do you want to continue, Y/N [Y]:','s'); 
if m_choice == 'n' || m_choice == 'N'
    break;
end;

% Example 2
% min x1 s.t. [x1;x2;x3] in Q(3)
% which can be written as 
%       min x1 
%       s.t. [0,0,0,1][x1;x2;x3;x4] = 1
%            [x1;x2;x3] in Q(3)
%                    x4 in R_+
clear blk At c;
blk{1,1}= 'q'; blk{1,2} = 3; At{1} = sparse(1,3); c{1} = [1;0;0];
blk{2,1} = 'l'; blk{2,2} = 1; At{2} = sparse(1); c{2} = 0;
b = 5;
[opt_obj, x_return, y_return, z_return, info] = hsd_lqeu_fast(blk, At, c,b);

% End of the infinite loop
display('End of all examples!');
break;
end
