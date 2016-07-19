function id_vec = get_id_soc(dim_list)
% Generate the identity element of a second order cone K = Q(q(1))×...×Q(q(n))
% Example: dim_list = [3;2;5] gives id_vec = [1;0;0;1;0;1;0;0;0;0]
id_vec = zeros(sum(dim_list),1);
for k = 1:length(dim_list)
    id_vec(sum(dim_list(1:k-1))+1) = 1;
end

