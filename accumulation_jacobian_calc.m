%% calculate the accumulation part of the Jacobian. 
function [accumulation_jacobian] = accumulation_jacobian_calc(P_vec, rock, fluid, grid, simulation)
delta_t = simulation.time_step;
sg = P_vec(2:2:size(P_vec));
phi = rock.porosity;

derivatives_blockcenter = derivatives_calc(P_vec, fluid, grid);
dbo_dpo = derivatives_blockcenter.dbo_dpo;
dbg_dpo = derivatives_blockcenter.dbg_dpo;
dRs_dpo = derivatives_blockcenter.dRs_dpo;

accumulation_jacobian = zeros(2*grid.blocknums);

for i = 1 : grid.blocknums 
    % accumulation term of the Jacobian on exist on the main diagonal blocks
    % for the gas equtaions
    accumulation_jacobian(2*i-1, 2*i-1) = (grid.Vij(i)/(5.615*delta_t))*((1-sg(i))...
        *(dRs_dpo(i)*phi*fluid.bo(i)+dbo_dpo(i)*fluid.Rs(i)*phi)+sg(i)*dbg_dpo(i)*phi);
    accumulation_jacobian(2*i-1, 2*i) = (grid.Vij(i)/(5.615*delta_t))*(-phi*fluid.Rs(i)*fluid.bo(i)+...
        phi*fluid.bg(i));
    % for the oil equtaions
    accumulation_jacobian(2*i, 2*i-1) = (grid.Vij(i)/(5.615*delta_t))*(1-sg(i))*(dbo_dpo(i)*phi);
    accumulation_jacobian(2*i, 2*i) = -(grid.Vij(i)/(5.615*delta_t))*(phi*fluid.bo(i));
end

accumulation_jacobian = sparse(accumulation_jacobian);


end