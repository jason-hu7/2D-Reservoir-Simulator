%% construct the derivatives matrix for Jacobian calculation
function [derivatives_matrix] = derivatives_matrix_calc(P_vec, grid, rock, fluid, diagonal)

pressure = P_vec(1:2:size(P_vec));
pressure_matrix = reshape(pressure, grid.Nx, grid.Ny)';
bo_matrix = reshape(fluid.bo, grid.Nx, grid.Ny)';
bg_matrix = reshape(fluid.bg, grid.Nx, grid.Ny)';
kro_matrix = reshape(fluid.k_ro, grid.Nx, grid.Ny)';
krg_matrix = reshape(fluid.k_rg, grid.Nx, grid.Ny)';
mu_g_matrix = reshape(fluid.gasViscosity, grid.Nx, grid.Ny)';

rock_geo_trans = rock_geo_transmissbility(grid, rock);
% derivatives of bo, bg, kro, mu_g, Rs with respect to primary variable at
% the center of each block
derivatives_blockcenter = derivatives_calc(P_vec, fluid, grid);
dbo_matrix = reshape(derivatives_blockcenter.dbo_dpo, grid.Nx, grid.Ny)';
dbg_matrix = reshape(derivatives_blockcenter.dbg_dpo, grid.Nx, grid.Ny)';
dkro_matrix = reshape(derivatives_blockcenter.dkro_dsg, grid.Nx, grid.Ny)';
dkrg_matrix = reshape(derivatives_blockcenter.dkrg_dsg, grid.Nx, grid.Ny)';
dmu_g_matrix = reshape(derivatives_blockcenter.dmu_g_dpo, grid.Nx, grid.Ny)';
dRs_matrix = reshape(derivatives_blockcenter.dRs_dpo, grid.Nx, grid.Ny)';

% 3rd dimension, 1 = left, 2 = right, 3 = above, 4 = below, 5 = center;
dgammao_dpo = zeros(grid.Ny, grid.Nx, 4);
dgammao_dsg = zeros(grid.Ny, grid.Nx, 4);
dgammag_dpo = zeros(grid.Ny, grid.Nx, 4);
dgammag_dsg = zeros(grid.Ny, grid.Nx, 4);
dRs_dpo = zeros(grid.Ny, grid.Nx, 4);

%% populate the matrix from col to col and then from row to row
% transmissbility from col to col
for col = 1 : grid.Nx
    if col == 1
        % derivatives of transmissbility to the left (1) is 0
        dgammao_dpo(:,col,1) = 0;
        dgammao_dsg(:,col,1) = 0;
        dgammag_dpo(:,col,1) = 0;
        dgammag_dsg(:,col,1) = 0;
        dRs_dpo(:,col,1) = 0;
        % calculate derivatives of transmissbility to the right
        Dp = zeros(grid.Ny,1); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(:,col) - pressure_matrix(:,col+1);
            dir_ind = 1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(:,col+1) - pressure_matrix(:,col);
            dir_ind = 0;
        end
        derivative_upwinding_col(Dp, col, dir_ind, 2);
        
    elseif col == grid.Nx
        % derivatives of transmissbility to the right (2) is 0
        dgammao_dpo(:,col,2) = 0;
        dgammao_dsg(:,col,2) = 0;
        dgammag_dpo(:,col,2) = 0;
        dgammag_dsg(:,col,2) = 0; 
        dRs_dpo(:,col,2) = 0;
        % calculate derivatives of transmissbility to the right
        Dp = zeros(grid.Ny,1); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(:,col) - pressure_matrix(:,col-1);
            dir_ind = -1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(:,col-1) - pressure_matrix(:,col);
            dir_ind = 0;
        end
        derivative_upwinding_col(Dp, col, dir_ind, 1);
        
    else
        % if the block is not on the edge then go both left and rght dir
        % calculate derivatives of transmissbility to the right
        Dp = zeros(grid.Ny,1); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(:,col) - pressure_matrix(:,col+1);
            dir_ind = 1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(:,col+1) - pressure_matrix(:,col);
            dir_ind = 0;
        end
        derivative_upwinding_col(Dp, col, dir_ind, 2);
        
        % now to the left direction
        Dp = zeros(grid.Ny,1); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(:,col) - pressure_matrix(:,col-1);
            dir_ind = -1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(:,col-1) - pressure_matrix(:,col);
            dir_ind = 0;
        end
        derivative_upwinding_col(Dp, col, dir_ind, 1);
    end
end

% transmissbility from row to row
for row = 1 : grid.Ny
    if row == 1
        % transmissbility to the up is (3) is 0
        dgammao_dpo(row,:,3) = 0;
        dgammao_dsg(row,:,3) = 0;
        dgammag_dpo(row,:,3) = 0;
        dgammag_dsg(row,:,3) = 0;
        dRs_dpo(row,:,3) = 0;
        % calculate derivatives of transmissbility to below
        Dp = zeros(1, grid.Nx); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(row,:) - pressure_matrix(row+1,:);
            dir_ind = 1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(row+1,:) - pressure_matrix(row,:);
            dir_ind = 0;
        end
        derivative_upwinding_row(Dp, row, dir_ind, 4);

    elseif row == grid.Ny
        % transmissbility to the down is (4) is 0
        dgammao_dpo(row,:,4) = 0;
        dgammao_dsg(row,:,4) = 0;
        dgammag_dpo(row,:,4) = 0;
        dgammag_dsg(row,:,4) = 0;
        dRs_dpo(row,:,4) = 0;
        % calculate derivatives of transmissbility to above (3)
        Dp = zeros(1, grid.Nx); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(row,:) - pressure_matrix(row-1,:);
            dir_ind = -1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(row-1,:) - pressure_matrix(row,:);
            dir_ind = 0;
        end
        derivative_upwinding_row(Dp, row, dir_ind, 3);
    else
         % calculate derivatives of transmissbility to above (3)
        Dp = zeros(1, grid.Nx); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(row,:) - pressure_matrix(row-1,:);
            dir_ind = -1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(row-1,:) - pressure_matrix(row,:);
            dir_ind = 0;
        end
        derivative_upwinding_row(Dp, row, dir_ind, 3);
        % calculate derivatives of transmissbility to below
        Dp = zeros(1, grid.Nx); dir_ind = 0;
        if strcmp(diagonal, 'non-main')
            Dp = pressure_matrix(row,:) - pressure_matrix(row+1,:);
            dir_ind = 1;
        elseif strcmp(diagonal, 'main')
            Dp = pressure_matrix(row+1,:) - pressure_matrix(row,:);
            dir_ind = 0;
        end
        derivative_upwinding_row(Dp, row, dir_ind, 4);
    end
end

derivatives_matrix.dgammag_dsg = dgammag_dsg;
derivatives_matrix.dgammag_dpo = dgammag_dpo;
derivatives_matrix.dgammao_dsg = dgammao_dsg;
derivatives_matrix.dgammao_dpo = dgammao_dpo;
derivatives_matrix.dRs_dpo = dRs_dpo;

    % calculate and upwind derivatives column by column
    function []= derivative_upwinding_col(Dp, col, dir_ind, direction)
        d_gammao_dpo_col = rock_geo_trans(:,col,direction).*kro_matrix(:,col+dir_ind).*dbo_matrix(:,col+dir_ind)./fluid.oilViscosity;
        d_gammao_dpo_col(Dp > 0) = 0;
        dgammao_dpo(:,col,direction) = d_gammao_dpo_col;
        
        d_gammao_dsg_col = rock_geo_trans(:,col,direction).*bo_matrix(:,col+dir_ind).*dkro_matrix(:,col+dir_ind)./fluid.oilViscosity;
        d_gammao_dsg_col(Dp > 0) = 0;
        dgammao_dsg(:,col,direction) = d_gammao_dsg_col;
        
        d_gammag_dpo_col = rock_geo_trans(:,col,direction).*(krg_matrix(:,col+dir_ind).*dbg_matrix(:,col+dir_ind)./mu_g_matrix(:,col+dir_ind) ...
            -krg_matrix(:,col+dir_ind).*bg_matrix(:,col+dir_ind).*dmu_g_matrix(:,col+dir_ind)./(mu_g_matrix(:,col+dir_ind).^2));
        d_gammag_dpo_col(Dp > 0) = 0;
        dgammag_dpo(:,col,direction) = d_gammag_dpo_col;
        
        d_gammag_dsg_col = rock_geo_trans(:,col,direction).*(bg_matrix(:,col+dir_ind).*dkrg_matrix(:,col+dir_ind)./mu_g_matrix(:,col+dir_ind));
        d_gammag_dsg_col(Dp > 0) = 0;
        dgammag_dsg(:,col,direction) = d_gammag_dsg_col;
        
        dRs_dpo_col = dRs_matrix(:,col+dir_ind);
        dRs_dpo_col(Dp > 0) = 0;
        dRs_dpo(:,col,direction) = dRs_dpo_col;
    end

    % calculate and upwind derivatives row by row
    function []= derivative_upwinding_row(Dp, row, dir_ind, direction)
        d_gammao_dpo_row = rock_geo_trans(row,:,direction).*kro_matrix(row+dir_ind,:).*dbo_matrix(row+dir_ind,:)./fluid.oilViscosity;
        d_gammao_dpo_row(Dp > 0) = 0;
        dgammao_dpo(row,:,direction) = d_gammao_dpo_row;
        
        d_gammao_dsg_row = rock_geo_trans(row,:,direction).*bo_matrix(row+dir_ind,:).*dkro_matrix(row+dir_ind,:)./fluid.oilViscosity;
        d_gammao_dsg_row(Dp > 0) = 0;
        dgammao_dsg(row,:,direction) = d_gammao_dsg_row;
        
        d_gammag_dpo_row = rock_geo_trans(row,:,direction).*(krg_matrix(row+dir_ind,:).*dbg_matrix(row+dir_ind,:)./mu_g_matrix(row+dir_ind,:) ...
            -krg_matrix(row+dir_ind,:).*bg_matrix(row+dir_ind,:).*dmu_g_matrix(row+dir_ind,:)./(mu_g_matrix(row+dir_ind,:).^2));
        d_gammag_dpo_row(Dp > 0) = 0;
        dgammag_dpo(row,:,direction) = d_gammag_dpo_row;
        
        d_gammag_dsg_row = rock_geo_trans(row,:,direction).*(bg_matrix(row+dir_ind,:).*dkrg_matrix(row+dir_ind,:)./mu_g_matrix(row+dir_ind,:));
        d_gammag_dsg_row(Dp > 0) = 0;
        dgammag_dsg(row,:,direction) = d_gammag_dsg_row;
        
        dRs_dpo_row = dRs_matrix(row+dir_ind,:);
        dRs_dpo_row(Dp > 0) = 0;
        dRs_dpo(row,:,direction) = dRs_dpo_row;
    end
    

end