%% calculate the fluid part of transmissbiulity
function [fluid_trans] = fluid_transmissbility(P_vec, fluid, grid)

% pressure = P(1:2:size(P));
% 3rd dimension, 1 = left, 2 = right, 3 = above, 4 = below;
fluid_trans_oil = zeros(grid.Ny, grid.Nx, 4);
fluid_trans_gas = zeros(grid.Ny, grid.Nx, 4);
fluid_trans_Rs = zeros(grid.Ny, grid.Nx, 4);

% transmissbility from col to col
for col = 1 : grid.Nx
    if col == 1
        % transmissbility to the left (1) is 0
        fluid_trans_oil(:,col,1) = 0;
        fluid_trans_gas(:,col,1) = 0;
        fluid_trans_Rs(:,col,1) = 0;
        right_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'right',col);
        fluid_trans_oil(:,col,2) = right_fluid_trans.oil;
        fluid_trans_gas(:,col,2) = right_fluid_trans.gas;
        fluid_trans_Rs(:,col,2) = right_fluid_trans.Rs;
    elseif col == grid.Nx
        % transmissbility to the right (2) is 0
        fluid_trans_oil(:,col,2) = 0;
        fluid_trans_gas(:,col,2) = 0;
        fluid_trans_Rs(:,col,2) = 0;
        left_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'left',col);
        fluid_trans_oil(:,col,1) = left_fluid_trans.oil;
        fluid_trans_gas(:,col,1) = left_fluid_trans.gas;
        fluid_trans_Rs(:,col,1) = left_fluid_trans.Rs;
    else
        % if the block is not on the edge then go both left and rght dir
        left_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'left',col);
        right_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'right',col);
        fluid_trans_oil(:,col,1) = left_fluid_trans.oil;
        fluid_trans_gas(:,col,1) = left_fluid_trans.gas;
        fluid_trans_Rs(:,col,1) = left_fluid_trans.Rs;
        fluid_trans_oil(:,col,2) = right_fluid_trans.oil;
        fluid_trans_gas(:,col,2) = right_fluid_trans.gas;
        fluid_trans_Rs(:,col,2) = right_fluid_trans.Rs;
    end
end

% transmissbility from row to row
for row = 1 : grid.Ny
    if row == 1
        % transmissbility to the up is (3) is 0
        fluid_trans_oil(row,:,3) = 0;
        fluid_trans_gas(row,:,3) = 0;
        fluid_trans_Rs(row,:,3) = 0;
        down_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'down',row);
        fluid_trans_oil(row,:,4) = down_fluid_trans.oil;
        fluid_trans_gas(row,:,4) = down_fluid_trans.gas;
        fluid_trans_Rs(row,:,4) = down_fluid_trans.Rs;
    elseif row == grid.Ny
        % transmissbility to the down is (4) is 0
        fluid_trans_oil(row,:,4) = 0;
        fluid_trans_gas(row,:,4) = 0;
        fluid_trans_Rs(row,:,4) = 0;
        up_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'up',row);
        fluid_trans_oil(row,:,3) = up_fluid_trans.oil;
        fluid_trans_gas(row,:,3) = up_fluid_trans.gas;
        fluid_trans_Rs(row,:,3) = up_fluid_trans.Rs;
    else
        up_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'up',row);
        down_fluid_trans = fluid_trans_upwind(P_vec,fluid,grid,'down',row);
        fluid_trans_oil(row,:,3) = up_fluid_trans.oil;
        fluid_trans_gas(row,:,3) = up_fluid_trans.gas;
        fluid_trans_Rs(row,:,3) = up_fluid_trans.Rs;
        fluid_trans_oil(row,:,4) = down_fluid_trans.oil;
        fluid_trans_gas(row,:,4) = down_fluid_trans.gas; 
        fluid_trans_Rs(row,:,4) = down_fluid_trans.Rs; 
    end
end

fluid_trans.oil = fluid_trans_oil;
fluid_trans.gas = fluid_trans_gas;
fluid_trans.Rs = fluid_trans_Rs;
end