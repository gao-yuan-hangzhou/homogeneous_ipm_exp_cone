# homogeneous_ipm_exp_cone


The main program is hsd_lqeu.m which takes in SDPT3-style cell array inputs. It solves problems of the following form:


min c'x

s.t. Ax = b, x(i) in K(i)

where x = [x(1); ..., x(N)], c = [c(1), ..., c(N)], A = [A(1), ..., A(N)]

and each K(j) can be one of the following:

  % (R_+)^nl = {x: x>=0}.

  % product of lorentz cones Q(n) = {x in R^n: x(1)>=||x(2:n)||}, where || || is the Euclidean 2-norm.

  % Product of the exponential cone K_exp = closure{(x1,x2,x3): x2>=0, x3>0, exp(x1/x3)<=x2/x3}.

  % R^nu for some integer nu.

How to call: 

[obj_val, x_return,y_return,z_return, info] = hsd_lqeu(blk, A_cell, c_cell, b, 1e-8, 500);

====== What are blk, A_cell, c_cell and b? ======

For j = 1,...,N

if K(j) is (R_+)^nl or R^nu, one has blk{j,1} = 'l' or 'u' respectively, blk{j,2} = dimension of x(j)

if K(j) is a product of Lorentz cones, one has blk{j,1} = 'q', blk{j,2} is, say, [5;3;2;5] (each element is the dimension of 
a single Lorenz cone)

if K(j) is a product of exponential cones, one has blk{j,1} = 'e', blk{j,2} = [3;3;...;3], length(blk{j,2}) = number of exponential cones in the product

Please refer to test_simple_examples, test_hsd_lqeu and test_hsd_lqeu_against_SDPT3 for more examples of usage.

Please email queries to gaoyuan@u.nus.edu or comment on the webpage. Thank you.

