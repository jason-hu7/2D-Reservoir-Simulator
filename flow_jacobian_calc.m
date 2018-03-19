%% construct the jacobian of the oil gas system for Newton's iterations
function [flow_jacobian] = flow_jacobian_calc(P_vec, rock, fluid, grid)

%% expand the derivatives calculated for transmissibilities and Rs for the non-main diagonal
derivatives_matrix = derivatives_matrix_calc(P_vec, grid, rock, fluid, 'non-main');
dgammao_dpo = derivatives_matrix.dgammao_dpo;
dgammao_dsg = derivatives_matrix.dgammao_dsg;
dgammag_dpo = derivatives_matrix.dgammag_dpo;
dgammag_dsg = derivatives_matrix.dgammag_dsg;
dRs_dpo = derivatives_matrix.dRs_dpo;

dgammao_dpo_right = dgammao_dpo(:,:,2)';
dgammao_dpo_left = dgammao_dpo(:,:,1)';
dgammao_dpo_up = dgammao_dpo(:,:,3)';
dgammao_dpo_down = dgammao_dpo(:,:,4)';

dgammao_dsg_right = dgammao_dsg(:,:,2)';
dgammao_dsg_left = dgammao_dsg(:,:,1)';
dgammao_dsg_up = dgammao_dsg(:,:,3)';
dgammao_dsg_down = dgammao_dsg(:,:,4)';

dgammag_dpo_right = dgammag_dpo(:,:,2)';
dgammag_dpo_left = dgammag_dpo(:,:,1)';
dgammag_dpo_up = dgammag_dpo(:,:,3)';
dgammag_dpo_down = dgammag_dpo(:,:,4)';

dgammag_dsg_right = dgammag_dsg(:,:,2)';
dgammag_dsg_left = dgammag_dsg(:,:,1)';
dgammag_dsg_up = dgammag_dsg(:,:,3)';
dgammag_dsg_down = dgammag_dsg(:,:,4)';

dRs_dpo_right = dRs_dpo(:,:,2)';
dRs_dpo_left = dRs_dpo(:,:,1)';
dRs_dpo_up = dRs_dpo(:,:,3)';
dRs_dpo_down = dRs_dpo(:,:,4)';
%-------------------------------------------------------------------------------------------%
%% expand the derivatives calculated for transmissibilities and Rs for the main diagonal
derivatives_matrix_main = derivatives_matrix_calc(P_vec, grid, rock, fluid, 'main');
dgammao_dpo_m = derivatives_matrix_main.dgammao_dpo;
dgammao_dsg_m = derivatives_matrix_main.dgammao_dsg;
dgammag_dpo_m = derivatives_matrix_main.dgammag_dpo;
dgammag_dsg_m = derivatives_matrix_main.dgammag_dsg;
dRs_dpo_m = derivatives_matrix_main.dRs_dpo;

dgammao_dpo_right_m = dgammao_dpo_m(:,:,2)';
dgammao_dpo_left_m = dgammao_dpo_m(:,:,1)';
dgammao_dpo_up_m = dgammao_dpo_m(:,:,3)';
dgammao_dpo_down_m = dgammao_dpo_m(:,:,4)';

dgammao_dsg_right_m = dgammao_dsg_m(:,:,2)';
dgammao_dsg_left_m = dgammao_dsg_m(:,:,1)';
dgammao_dsg_up_m = dgammao_dsg_m(:,:,3)';
dgammao_dsg_down_m = dgammao_dsg_m(:,:,4)';

dgammag_dpo_right_m = dgammag_dpo_m(:,:,2)';
dgammag_dpo_left_m = dgammag_dpo_m(:,:,1)';
dgammag_dpo_up_m = dgammag_dpo_m(:,:,3)';
dgammag_dpo_down_m = dgammag_dpo_m(:,:,4)';

dgammag_dsg_right_m = dgammag_dsg_m(:,:,2)';
dgammag_dsg_left_m = dgammag_dsg_m(:,:,1)';
dgammag_dsg_up_m = dgammag_dsg_m(:,:,3)';
dgammag_dsg_down_m = dgammag_dsg_m(:,:,4)';

dRs_dpo_right_m = dRs_dpo_m(:,:,2)';
dRs_dpo_left_m = dRs_dpo_m(:,:,1)';
dRs_dpo_up_m = dRs_dpo_m(:,:,3)';
dRs_dpo_down_m = dRs_dpo_m(:,:,4)';
%--------------------------------------------------------------------------------%
%% construct the Dp matrix for later calculation
Dp_matrix = Dp_matrix_construct(P_vec, grid);
Dp_right = Dp_matrix(:,:,2)';
Dp_left = Dp_matrix(:,:,1)';
Dp_up = Dp_matrix(:,:,3)';
Dp_down = Dp_matrix(:,:,4)';

%% import the rock and geometry transmissibility and upwinded fluid transmissbility for calculation
fluid_trans = fluid_transmissbility(P_vec, fluid, grid);
rock_geo_trans = rock_geo_transmissbility(grid, rock);

gammao_right = rock_geo_trans(:,:,2)' .* fluid_trans.oil(:,:,2)';
gammao_left = rock_geo_trans(:,:,1)' .* fluid_trans.oil(:,:,1)';
gammao_up = rock_geo_trans(:,:,3)' .* fluid_trans.oil(:,:,3)';
gammao_down = rock_geo_trans(:,:,4)' .* fluid_trans.oil(:,:,4)';

gammag_right = rock_geo_trans(:,:,2)' .* fluid_trans.gas(:,:,2)';
gammag_left = rock_geo_trans(:,:,1)' .* fluid_trans.gas(:,:,1)';
gammag_up = rock_geo_trans(:,:,3)' .* fluid_trans.gas(:,:,3)';
gammag_down = rock_geo_trans(:,:,4)' .* fluid_trans.gas(:,:,4)';

Rs_right = fluid_trans.Rs(:,:,2)';
Rs_left = fluid_trans.Rs(:,:,1)';
Rs_up = fluid_trans.Rs(:,:,3)';
Rs_down = fluid_trans.Rs(:,:,4)';
%-------------------------------------------------------------------------------%
%% calculate the Jacobian
flow_jacobian = zeros(2*grid.blocknums);

for i = 1 : grid.blocknums
    % for the i+1 diagonal
    if i ~= grid.blocknums
        % for the gas equations
        flow_jacobian(2*i-1, 2*i+1) = gammag_right(i)+Rs_right(i).*gammao_right(i)+...
            (dgammag_dpo_right(i)+Rs_right(i).*dgammao_dpo_right(i)+dRs_dpo_right(i).*gammao_right(i)).*Dp_right(i);
        flow_jacobian(2*i-1, 2*i+2) = (dgammag_dsg_right(i)+Rs_right(i).*dgammao_dsg_right(i)).*Dp_right(i);
        % for the oil equations
        flow_jacobian(2*i, 2*i+1) = gammao_right(i) + dgammao_dpo_right(i).*Dp_right(i);
        flow_jacobian(2*i, 2*i+2) = dgammao_dsg_right(i).*Dp_right(i);
    end
    % for the i-1 diagonal
    if i ~= 1
        % for the gas equations
        flow_jacobian(2*i-1, 2*i-3) = gammag_left(i)+Rs_left(i).*gammao_left(i)+...
            (dgammag_dpo_left(i)+Rs_left(i).*dgammao_dpo_left(i)+dRs_dpo_left(i).*gammao_left(i)).*Dp_left(i);
        flow_jacobian(2*i-1, 2*i-2)= (dgammag_dsg_left(i)+Rs_left(i).*dgammao_dsg_left(i)).*Dp_left(i);
        % for the oil equations
        flow_jacobian(2*i, 2*i-3) = gammao_left(i) + dgammao_dpo_left(i).*Dp_left(i);
        flow_jacobian(2*i, 2*i-2) = dgammao_dsg_left(i).*Dp_left(i);
    end
    % for the i+Nx diagonal
    if i <= grid.blocknums - grid.Nx
        % for the gas equations
        flow_jacobian(2*i-1, 2*(i+grid.Nx)-1) = gammag_down(i)+Rs_down(i).*gammao_down(i)+...
            (dgammag_dpo_down(i)+Rs_down(i).*dgammao_dpo_down(i)+dRs_dpo_down(i).*gammao_down(i)).*Dp_down(i);
        flow_jacobian(2*i-1, 2*(i+grid.Nx)) = (dgammag_dsg_down(i)+Rs_down(i).*dgammao_dsg_down(i)).*Dp_down(i);
        % for the oil equations
        flow_jacobian(2*i, 2*(i+grid.Nx)-1) = gammao_down(i) + dgammao_dpo_down(i).*Dp_down(i);
        flow_jacobian(2*i, 2*(i+grid.Nx)) = dgammao_dsg_down(i).*Dp_down(i);
    end
    % for the i-Nx diagonal
    if i  > grid.Nx
        % for the gas equtaions
        flow_jacobian(2*i-1, 2*(i-grid.Nx)-1) = gammag_up(i)+Rs_up(i).*gammao_up(i)+...
            (dgammag_dpo_up(i)+Rs_up(i).*dgammao_dpo_up(i)+dRs_dpo_up(i).*gammao_up(i)).*Dp_up(i);
        flow_jacobian(2*i-1, 2*(i-grid.Nx)) =(dgammag_dsg_up(i)+Rs_up(i).*dgammao_dsg_up(i)).*Dp_up(i);
        % for the oil equations
        flow_jacobian(2*i, 2*(i-grid.Nx)-1) = gammao_up(i) + dgammao_dpo_up(i).*Dp_up(i);
        flow_jacobian(2*i, 2*(i-grid.Nx)) = dgammao_dsg_up(i).*Dp_up(i);
    end
    % for the main diagonal block
    % for the gas equtaions
    flow_jacobian(2*i-1, 2*i-1) = (dgammag_dpo_up_m(i)+dRs_dpo_up_m(i).*gammao_up(i)+dgammao_dpo_up_m(i).*Rs_up(i)).*Dp_up(i)...
        +(dgammag_dpo_left_m(i)+dRs_dpo_left_m(i).*gammao_left(i)+dgammao_dpo_left_m(i).*Rs_left(i)).*Dp_left(i)...
        +(dgammag_dpo_down_m(i)+dRs_dpo_down_m(i).*gammao_down(i)+dgammao_dpo_down_m(i).*Rs_down(i)).*Dp_down(i)...
        +(dgammag_dpo_right_m(i)+dRs_dpo_right_m(i).*gammao_right(i)+dgammao_dpo_right_m(i).*Rs_right(i)).*Dp_right(i)...
        -(gammag_up(i)+gammag_down(i)+gammag_left(i)+gammag_right(i))...
        -(Rs_up(i).*gammao_up(i)+Rs_left(i).*gammao_left(i)+Rs_down(i).*gammao_down(i)+Rs_right(i).*gammao_right(i));
    flow_jacobian(2*i-1, 2*i) = (dgammag_dsg_up_m(i)+dgammao_dsg_up_m(i).*Rs_up(i)).*Dp_up(i)...
        +(dgammag_dsg_left_m(i)+dgammao_dsg_left_m(i).*Rs_left(i)).*Dp_left(i)...
        +(dgammag_dsg_down_m(i)+dgammao_dsg_down_m(i).*Rs_down(i)).*Dp_down(i)...
        +(dgammag_dsg_right_m(i)+dgammao_dsg_right_m(i).*Rs_right(i)).*Dp_right(i);
    % for the oil equtaions
    flow_jacobian(2*i, 2*i-1) = dgammao_dpo_up_m(i).*Dp_up(i)+dgammao_dpo_down_m(i).*Dp_down(i)...
        +dgammao_dpo_left_m(i).*Dp_left(i)+dgammao_dpo_right_m(i).*Dp_right(i)...
        -(gammao_up(i)+gammao_down(i)+gammao_left(i)+gammao_right(i));
    flow_jacobian(2*i, 2*i) = dgammao_dsg_up_m(i).*Dp_up(i)+dgammao_dsg_left_m(i).*Dp_left(i)...
        +dgammao_dsg_right_m(i).*Dp_right(i)+dgammao_dsg_down_m(i).*Dp_down(i);
end

flow_jacobian = sparse(flow_jacobian);

end