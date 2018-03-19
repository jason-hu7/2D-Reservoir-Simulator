%% use upwinding to treat fluid part of transmissbility between blocks
function [fluid_trans] = fluid_trans_upwind(P_vec, fluid, grid, direction, col_row_num) 

pressure = P_vec(1:2:size(P_vec));
pressure_matrix = reshape(pressure, grid.Nx, grid.Ny)';
bo_matrix = reshape(fluid.bo, grid.Nx, grid.Ny)';
bg_matrix = reshape(fluid.bg, grid.Nx, grid.Ny)';
kro_matrix = reshape(fluid.k_ro, grid.Nx, grid.Ny)';
krg_matrix = reshape(fluid.k_rg, grid.Nx, grid.Ny)';
mu_gas_matrix = reshape(fluid.gasViscosity, grid.Nx, grid.Ny)';
Rs_matrix = reshape(fluid.Rs, grid.Nx, grid.Ny)';

fluid_trans_calc = @(kr, b, mu) (kr .* b ./ mu);

% for left direction
if strcmp('left',direction)
    ind = -1;
    Dp = pressure_matrix(:,col_row_num) - pressure_matrix(:,col_row_num-1);
    fluid_trans = x_dir_trans_calc(Dp,fluid,grid,col_row_num,ind);
    
elseif strcmp('right',direction)
    ind = 1;
    Dp = pressure_matrix(:,col_row_num) - pressure_matrix(:,col_row_num+1);
    fluid_trans = x_dir_trans_calc(Dp,fluid,grid,col_row_num,ind);
    
elseif strcmp('up',direction)
    ind = -1;
    Dp = pressure_matrix(col_row_num,:) - pressure_matrix(col_row_num-1,:);
    fluid_trans = y_dir_trans_calc(Dp,fluid,grid,col_row_num,ind);
    
elseif strcmp('down',direction)
    ind = 1;
    Dp = pressure_matrix(col_row_num,:) - pressure_matrix(col_row_num+1,:);
    fluid_trans = y_dir_trans_calc(Dp,fluid,grid,col_row_num,ind);   
end

    function [trans] = x_dir_trans_calc(Dp,fluid,grid,col_num,ind)    
        trans.oil = zeros(grid.Ny, 1);
        trans.gas = zeros(grid.Ny, 1);
        trans.Rs = zeros(grid.Ny, 1);
        
        kro_current_col = kro_matrix(:,col_num);
        kro_neighbor_col = kro_matrix(:,col_num+ind);   
        bo_current_col = bo_matrix(:,col_num);
        bo_neighbor_col = bo_matrix(:,col_num+ind);
        neighbor_trans_oil(Dp>=0) = fluid_trans_calc(kro_current_col(Dp>=0), bo_current_col(Dp>=0), fluid.oilViscosity);
        neighbor_trans_oil(Dp<0) = fluid_trans_calc(kro_neighbor_col(Dp<0), bo_neighbor_col(Dp<0), fluid.oilViscosity);
        trans.oil = neighbor_trans_oil;

        krg_current_col = krg_matrix(:,col_num);
        krg_neighbor_col = krg_matrix(:,col_num+ind);
        bg_current_col = bg_matrix(:,col_num);
        bg_neighbor_col = bg_matrix(:,col_num+ind);
        mu_current_col = mu_gas_matrix(:,col_num);
        mu_neighbor_col = mu_gas_matrix(:,col_num+ind);
        neighbor_trans_gas(Dp>=0) = fluid_trans_calc(krg_current_col(Dp>=0), bg_current_col(Dp>=0), mu_current_col(Dp>=0));
        neighbor_trans_gas(Dp<0) = fluid_trans_calc(krg_neighbor_col(Dp<0), bg_neighbor_col(Dp<0), mu_neighbor_col(Dp<0));
        trans.gas = neighbor_trans_gas;
        
        Rs_current_col = Rs_matrix(:,col_num);
        Rs_neighbor_col = Rs_matrix(:,col_num+ind);
        neighbor_Rs(Dp>=0) = Rs_current_col(Dp>=0);
        neighbor_Rs(Dp<0) = Rs_neighbor_col(Dp<0);
        trans.Rs = neighbor_Rs;
    end

    function trans = y_dir_trans_calc(Dp,fluid,grid,row_num,ind)
        trans.oil = zeros(1,grid.Nx);
        trans.gas = zeros(1,grid.Nx);
        
        kro_current_row = kro_matrix(row_num,:);
        kro_neighbor_row = kro_matrix(row_num+ind,:);   
        bo_current_row = bo_matrix(row_num,:);
        bo_neighbor_row = bo_matrix(row_num+ind,:);
        neighbor_trans_oil(Dp>=0) = fluid_trans_calc(kro_current_row(Dp>=0), bo_current_row(Dp>=0), fluid.oilViscosity);
        neighbor_trans_oil(Dp<0) = fluid_trans_calc(kro_neighbor_row(Dp<0), bo_neighbor_row(Dp<0), fluid.oilViscosity);
        trans.oil = neighbor_trans_oil;
        
        krg_current_row = krg_matrix(row_num,:);
        krg_neighbor_row = krg_matrix(row_num+ind,:);
        bg_current_row = bg_matrix(row_num,:);
        bg_neighbor_row = bg_matrix(row_num+ind,:);
        mu_current_row = mu_gas_matrix(row_num,:);
        mu_neighbor_row = mu_gas_matrix(row_num+ind,:);
        neighbor_trans_gas(Dp>=0) = fluid_trans_calc(krg_current_row(Dp>=0), bg_current_row(Dp>=0), mu_current_row(Dp>=0));
        neighbor_trans_gas(Dp<0) = fluid_trans_calc(krg_neighbor_row(Dp<0), bg_neighbor_row(Dp<0), mu_neighbor_row(Dp<0));
        trans.gas = neighbor_trans_gas;
        
        Rs_current_row = Rs_matrix(row_num,:);
        Rs_neighbor_row = Rs_matrix(row_num+ind,:);
        neighbor_Rs(Dp>=0) = Rs_current_row(Dp>=0);
        neighbor_Rs(Dp<0) = Rs_neighbor_row(Dp<0);
        trans.Rs = neighbor_Rs;
    end

end