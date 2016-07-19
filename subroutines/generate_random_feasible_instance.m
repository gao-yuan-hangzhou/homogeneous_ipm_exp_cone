function [A_cell, c_cell, b] = generate_random_feasible_instance(blk, m)
% This function takes in the dimension of the blocks of type 
% 'l', 'q', 'e' and 'u' and
% returns a (stictly) primal and dual feasible instance
% in terms of A_cell (cell array), c_cell (cell array), b (vector)

% Random sparse density function
density = @() 0.15+0.05*rand()-0.05*rand();

% Get the total number of blocks
N_block = size(blk,1);

% Get the total number of linear constraints or randomly set m < total_dimension_of_x
total_dim = 0;
for k = 1:N_block
    total_dim = total_dim + sum(blk{k,2});
end

if nargin < 2
    m = randi(total_dim);
end

% Generate A_cell
for k = 1:N_block
    A_cell{k} = sprandn(m, sum(blk{k,2}), density());
end

% Randomly generate x = [x_cone; x_unrestricted] and z = [z_cone; z_null]
% such that x_cone is in int(K), z_cone is in int(K^*), z_null = 0
for k = 1:N_block
    if blk{k,1} == 'l'
        x_cell{k} = exp(randn(blk{k,2},1)) .* exp(randn(blk{k,2},1)); 
        z_cell{k} = exp(randn(blk{k,2},1)) .* exp(randn(blk{k,2},1));
    elseif blk{k,1} == 'q'
        dim_q_array = blk{k,2};
        x_cell{k} = zeros(sum(dim_q_array),1);
        z_cell{k} = zeros(sum(dim_q_array),1);
        for j = 1:length(dim_q_array)
            curr_dim = dim_q_array(j);
            curr_begin = sum(dim_q_array(1:j-1))+1;
            curr_end = sum(dim_q_array(1:j));
            x_cell{k}(curr_begin+1:curr_end) = exp(randn()) * randn(curr_dim-1,1);
            x_cell{k}(curr_begin) = exp(randn()) + sqrt(x_cell{k}(curr_begin+1:curr_end)'*x_cell{k}(curr_begin+1:curr_end));
            z_cell{k}(curr_begin+1:curr_end) = exp(randn()) * randn(curr_dim-1,1);
            z_cell{k}(curr_begin) = exp(randn()) + sqrt(z_cell{k}(curr_begin+1:curr_end)'*z_cell{k}(curr_begin+1:curr_end));
        end
    elseif blk{k,1} == 'e'
        Ne = length(blk{k,2});
        exp_cone_center = [-1.0151; 1.2590; 0.5560];
        x_cell{k} = zeros(3*Ne,1);
        z_cell{k} = zeros(3*Ne,1);
        for j = 1:Ne
            x_cell{k}(3*j-2:3*j) = exp(randn())*exp_cone_center;
            z_cell{k}(3*j-2:3*j) = exp(randn())*exp_cone_center;
        end
    elseif blk{k,1} == 'u'
        x_cell{k} = exp(randn()) * randn(blk{k,2},1);
        z_cell{k} = zeros(blk{k,2},1);
    end
end

% Calculate b = Ax = A(1)x(1)+...+A(N)x(N) and c = A'y + z
b = zeros(m,1); y = exp(randn()) * randn(m,1);
for k = 1:N_block
    b = b + A_cell{k}*x_cell{k};
    c_cell{k} = A_cell{k}'*y + z_cell{k};
end

% Function ends
end

