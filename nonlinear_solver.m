%% use Newton's method to solve the simulation
function [P_n_k, simulation, well] = nonlinear_solver(P_n, fluid, rock, grid, rock_geo_trans, simulation, well)
k_max = 10;
k = 0;
P_n_k = P_n;

% calculate the fluid properties based on the P vector
while k <= k_max
    k = k + 1;
    fluid = fluid_property_calc(P_n_k, fluid, grid);

    % calculate the residual
    fluid_T_k = fluid_transmissbility(P_n_k, fluid, grid);
    T_matrix_k = transmissbility_mat_construct(grid, fluid_T_k, rock_geo_trans);
    A_k = accumulation_construct(P_n_k, P_n, rock, fluid, grid, simulation);
    [W_k,well] = well_construct(P_n_k, well, grid, fluid);
    R_k = residual_calc(T_matrix_k, P_n_k, A_k, W_k);
    
    % calculate the Jacobian
    J_f_k = flow_jacobian_calc(P_n_k, rock, fluid, grid);
    J_a_k = accumulation_jacobian_calc(P_n_k, rock, fluid, grid, simulation);
    J_w_k = well_jacobian_calc(P_n_k, well, fluid, grid);
    J_matrix_k = jacobian_matrix_construct(J_f_k, J_a_k, J_w_k);
    
    delta_k = - J_matrix_k \ R_k;
    P_n_k1 = P_n_k + delta_k;
    % fix negative sat and sat over 1
    for i = 2 : 2 : size(P_n_k1)
        if P_n_k1(i) < 0
            P_n_k1(i) = 0;
        elseif P_n_k1(i) > 1
            P_n_k1(i) = 1;
        end
    end
    
    % get new residual after update
    fluid = fluid_property_calc(P_n_k1, fluid, grid);
    fluid_T_k_1 = fluid_transmissbility(P_n_k1, fluid, grid);
    T_matrix_k_1 = transmissbility_mat_construct(grid, fluid_T_k_1, rock_geo_trans);
    A_k_1 = accumulation_construct(P_n_k1, P_n, rock, fluid, grid, simulation);
    [W_k_1,well] = well_construct(P_n_k1, well, grid, fluid);
    R_k_1 = residual_calc(T_matrix_k_1, P_n_k1, A_k_1, W_k_1);
    
    if convergence_check(R_k_1, P_n_k1, P_n_k, grid, rock, fluid, simulation)
        fprintf('converged! \n');
        break;
    end
    
    % cut time by half if doens't converge
    P_n_k = P_n_k1;
    if k > k_max
        fprintf('max newton loop iteration exceeded, cut time step by half and restart! \n')
        simulation.time_step = simulation.time_step / 2;
        k = 0;
        P_n_k = P_n; 
    end
    
end

simulation.newton_num = k;

end