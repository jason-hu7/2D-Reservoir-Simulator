%% put the jacobian together, will add wells term later
function [jacobian_matrix] = jacobian_matrix_construct(J_f, J_a, J_w)

    jacobian_matrix = J_f - J_a + J_w;

end