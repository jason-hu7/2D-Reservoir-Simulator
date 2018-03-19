%% calculate the rock and geometric part of transmissbility
function [rock_geo_trans] = rock_geo_transmissbility(grid, rock)

alpha = 0.001127;
% kx_matrix = reshape(rock.kx, grid.Nx, grid.Ny)';
% ky_matrix = reshape(rock.ky, grid.Nx, grid.Ny)';
kx_matrix = rock.kx;
ky_matrix = rock.ky;
% 3rd dimension, 1 = left, 2 = right, 3 = above, 4 = below;
rock_geo_trans = zeros(grid.Ny, grid.Nx, 4);

% calculate delta x for the first row of grid
delta_x_1stRow = (grid.Lx / grid.Nx) * ones(1, grid.Nx);
delta_x = repmat(delta_x_1stRow, grid.Ny, 1);

delta_y_1stCol = (grid.Ly / grid.Ny) * ones(grid.Ny, 1);
delta_y = repmat(delta_y_1stCol, 1, grid.Nx);

delta_z = (grid.Lz / grid.Nz) * ones(grid.Ny, grid.Nx);


%% rock and geo transmissbilit between one column to another column
k_harmonic_avg = @(l_1,l_2,k_1,k_2) (l_1 + l_2) ./ (l_1 ./ k_1 + l_2 ./ k_2);

for col = 1 : grid.Nx
    if col == 1
        rock_geo_trans(:,col,1) = 0;
        kx_avg_right = k_harmonic_avg(delta_x(:,col), delta_x(:,col+1), kx_matrix(:,col), kx_matrix(:,col+1));
        rock_geo_trans(:,col,2) = alpha * kx_avg_right .* delta_y(:,col) .* delta_z(:,col) ./ delta_x(:,col);
    elseif col == grid.Nx
        rock_geo_trans(:,col,2) = 0;
        kx_avg_left = k_harmonic_avg(delta_x(:,col), delta_x(:,col-1), kx_matrix(:,col), kx_matrix(:,col-1));
        rock_geo_trans(:,col,1) = alpha * kx_avg_left .* delta_y(:,col) .* delta_z(:,col) ./ delta_x(:,col);
    else
        % trans to left and right are all calculated
        kx_avg_left = k_harmonic_avg(delta_x(:,col), delta_x(:,col-1), kx_matrix(:,col), kx_matrix(:,col-1));
        kx_avg_right = k_harmonic_avg(delta_x(:,col), delta_x(:,col+1), kx_matrix(:,col), kx_matrix(:,col+1));
        rock_geo_trans(:,col,1) = alpha * kx_avg_left .* delta_y(:,col) .* delta_z(:,col) ./ delta_x(:,col);
        rock_geo_trans(:,col,2) = alpha * kx_avg_right .* delta_y(:,col) .* delta_z(:,col) ./ delta_x(:,col);
    end
end

%% rock and geo transmissbility between one row to another row
for row = 1 : grid.Ny
    if row == 1
        rock_geo_trans(row,:,3) = 0;
        ky_avg_down = k_harmonic_avg(delta_y(row,:), delta_y(row+1,:), ky_matrix(row,:), ky_matrix(row+1,:));
        rock_geo_trans(row,:,4) = alpha * ky_avg_down .* delta_x(row,:) .* delta_z(row,:) ./ delta_y(row,:);
    elseif row == grid.Ny
        rock_geo_trans(row,:,4) = 0;
        ky_avg_up = k_harmonic_avg(delta_y(row,:), delta_y(row-1,:), ky_matrix(row,:), ky_matrix(row-1,:));
        rock_geo_trans(row,:,3) = alpha * ky_avg_up .* delta_x(row,:) .* delta_z(row,:) ./ delta_y(row,:);
    else
        ky_avg_up = k_harmonic_avg(delta_y(row,:), delta_y(row-1,:), ky_matrix(row,:), ky_matrix(row-1,:));
        rock_geo_trans(row,:,3) = alpha * ky_avg_up .* delta_x(row,:) .* delta_z(row,:) ./ delta_y(row,:);
        ky_avg_down = k_harmonic_avg(delta_y(row,:), delta_y(row+1,:), ky_matrix(row,:), ky_matrix(row+1,:));
        rock_geo_trans(row,:,4) = alpha * ky_avg_down .* delta_x(row,:) .* delta_z(row,:) ./ delta_y(row,:);
    end
end

end