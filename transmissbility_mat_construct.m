%% construct the transmissibility matrix, this matrix will be a pentadiagonal matrix
function [transm_matrix] = transmissbility_mat_construct(grid, fluid_trans, rock_geo_trans)

gamma_ox_right = rock_geo_trans(:,:,2)' .* fluid_trans.oil(:,:,2)';
gamma_ox_left = rock_geo_trans(:,:,1)' .* fluid_trans.oil(:,:,1)';
gamma_oy_up = rock_geo_trans(:,:,3)' .* fluid_trans.oil(:,:,3)';
gamma_oy_down = rock_geo_trans(:,:,4)' .* fluid_trans.oil(:,:,4)';

gamma_gx_right = rock_geo_trans(:,:,2)' .* fluid_trans.gas(:,:,2)';
gamma_gx_left = rock_geo_trans(:,:,1)' .* fluid_trans.gas(:,:,1)';
gamma_gy_up = rock_geo_trans(:,:,3)' .* fluid_trans.gas(:,:,3)';
gamma_gy_down = rock_geo_trans(:,:,4)' .* fluid_trans.gas(:,:,4)';

Rs_gamma_ox_right = fluid_trans.Rs(:,:,2)' .* gamma_ox_right;
Rs_gamma_ox_left = fluid_trans.Rs(:,:,1)' .* gamma_ox_left;
Rs_gamma_oy_up = fluid_trans.Rs(:,:,3)' .* gamma_oy_up;
Rs_gamma_oy_down = fluid_trans.Rs(:,:,4)' .* gamma_oy_down;

transm_matrix = zeros(2 * grid.blocknums);

for i = 1 : grid.blocknums
    % for the i+1 diagonal
    if i ~= grid.blocknums
        transm_matrix(2*i-1, 2*i+1) = Rs_gamma_ox_right(i) + gamma_gx_right(i);
        transm_matrix(2*i, 2*i+1) = gamma_ox_right(i);
    end
    % for the i-1 diagonal
    if i ~= 1
        transm_matrix(2*i-1, 2*i-3) = Rs_gamma_ox_left(i) + gamma_gx_left(i);
        transm_matrix(2*i, 2*i-3) = gamma_ox_left(i);
    end
    % for the i+Nx diagonal
    if i <= (grid.blocknums - grid.Nx)
        transm_matrix(2*i-1, 2*(i+grid.Nx)-1) = Rs_gamma_oy_down(i) + gamma_gy_down(i);
        transm_matrix(2*i, 2*(i+grid.Nx)-1) = gamma_oy_down(i);
    end
    % for the i-Nx diagonal
    if i  > grid.Nx
        transm_matrix(2*i-1, 2*(i-grid.Nx)-1) = Rs_gamma_oy_up(i) + gamma_gy_up(i);
        transm_matrix(2*i, 2*(i-grid.Nx)-1) = gamma_oy_up(i);
    end
end

% sum up each row to get the element located at i,i
for i = 1 : grid.blocknums
    transm_matrix(2*i-1, 2*i-1) = -sum(transm_matrix(2*i-1,:));
    transm_matrix(2*i, 2*i-1) = -sum(transm_matrix(2*i,:));
end

transm_matrix = sparse(transm_matrix);

end