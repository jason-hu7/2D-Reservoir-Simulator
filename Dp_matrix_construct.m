%% calculate the pressure difference detween cell and its neighbors and store it in matrix form
function [Dp_matrix] = Dp_matrix_construct(P_vec, grid)

pressure = P_vec(1:2:size(P_vec));
pressure_matrix = reshape(pressure, grid.Nx, grid.Ny)';
Dp_matrix = zeros(grid.Ny, grid.Nx, 4);

%% calculate pressure difference from one column to another
for col = 1 : grid.Nx
    if col == 1
        Dp_matrix(:,col,1) = 0;
        Dp_matrix(:,col,2) = pressure_matrix(:,col+1)-pressure_matrix(:,col);
    elseif col == grid.Nx
        Dp_matrix(:,col,2) = 0;
        Dp_matrix(:,col,1) = pressure_matrix(:,col-1)-pressure_matrix(:,col);
    else
        Dp_matrix(:,col,1) = pressure_matrix(:,col-1)-pressure_matrix(:,col);
        Dp_matrix(:,col,2) = pressure_matrix(:,col+1)-pressure_matrix(:,col);
    end
end

%% calculate pressure difference from one row to another
for row = 1 : grid.Ny
    if row == 1
        Dp_matrix(row,:,3) = 0;
        Dp_matrix(row,:,4) = pressure_matrix(row+1,:)-pressure_matrix(row,:);
    elseif row == grid.Ny
        Dp_matrix(row,:,4) = 0;
        Dp_matrix(row,:,3) = pressure_matrix(row-1,:)-pressure_matrix(row,:);
    else
        Dp_matrix(row,:,3) = pressure_matrix(row-1,:)-pressure_matrix(row,:);
        Dp_matrix(row,:,4) = pressure_matrix(row+1,:)-pressure_matrix(row,:);
    end
end

end