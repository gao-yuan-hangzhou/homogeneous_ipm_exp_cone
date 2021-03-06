# homogeneous_ipm_exp_cone

Please email me at gaoyuan@u.nus.edu for the latest codes. 

Presentation slides are available at https://goo.gl/vTK0WA, which include several important examples involving exponential cone constraints and a brief description of the algorithm implemented. In short, it solves for the search directions though solving the Schur complement equations of the original linear systems using the preconditioned BiCGSTAB. The preconditioner is formed using the Sherman-Morrison-Woodburry formua applied on the Schur complement matrix, which is the sum of its sparse positive semidefinite part plus a low-rank perturbation term. This allows the solver to make use of the highly efficient and robust Cholesky factorization subroutine (on the sparse PSD part of the Schur complement matrix).

The main program is [hsd_lqeu.p] which takes in SDPT3-style (http://www.optimization-online.org/DB_FILE/2010/06/2654.pdf) cell array inputs. It solves problems coded in the following format:

min c{1}' * x{1} + ... + c{N}' * x{N}

s.t. A{1} * x{1}+...+A{n} * x{n} = b, x{j} ∈ K{j}, i = 1, ..., N

and each K{j}, j = 1, ..., N can be one of the following:

  Linear cone: blk{j,1} = 'l' means x{j} ∈ nonnegative orthant of dimension blk{j,2} or sum(blk{j,2})
  
  Second-order cone: blk{j,1} = 'q' means x{j} ∈ product of second-order cones. In this case, blk{j,2} = [q(1); ...; q(n)] gives the dimensions of the individual second-order cones. A second order cone of dimension p is defined as Q(p)= {x is p-dimensional: x(1)>=||x(2:n)||}, where ||·|| is the usual 2-norm.

  Exponential cone: blk{j,1} = 'e' means x{j} ∈ product of the exponential cone Kexp = closure{(x1,x2,x3): x2>=0, x3>0, exp(x1/x3)<=x2/x3}. In this case, there are length(blk{j,2}) exponential cones, where blk{j,2} = [3;...;3].

  Unrestricted space: blk{j,1} = 'u' means x{j} has dimension blk{j,2} and is unrestricted.


## ======== INPUT FORMAT ======== ##

[x_re,y_re,z_re, info_re] = hsd_lueq_fast(blk, A, c, b, input_options)

blk is the (2-by-N) cell array storing the dimensions of individual x{j}, as described above.

A is the cell array of coefficient matrices consisting of A{1}, ... , A{N}.

c is the cell array of cost vectors consisting of c{1}, ..., c{N}.

b is the right hand side vector in A{1} * x{1}+...+A{n} * x{n} = b.

You may wish to set one or more fields in the structure input_options as follows (note that everything below is OPTIONAL)

input_options.rel_eps: a number specifies the desired relative accuracy. If not specified, the default value is 1e-8.

input_options.max_iter_count: the maximum number of iteartions allowed. If not specified, the default value is 500.

input_options.initial_x: a cell array specifying x as part of the initial iterate.

input_options.initial_y: a vector specifying y as part of the initial iterate. Note that its dimension is m=dim(b).

input_options.initial_z: a cell array specifying z as part of the initial iterate.

## ======== References ======== ##

Please refer to the references in the presentation slides.
