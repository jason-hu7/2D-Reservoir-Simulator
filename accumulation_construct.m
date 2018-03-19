%% to construct the accumulation matrix D
function [accumulation] = accumulation_construct(P_n_1, P_n, rock, fluid, grid, simulation)
% oil pressure and gas saturation at time step n
delta_t = simulation.time_step;
po_n = P_n(1:2:size(P_n));
sg_n = P_n(2:2:size(P_n));
% oil pressure and gas saturation at time step n + 1
po_n_1 = P_n_1(1:2:size(P_n_1));
sg_n_1 = P_n_1(2:2:size(P_n_1));
% initiate the accumulation matrix
accum_vec = zeros(2 * grid.blocknums, 1);

FVF_n = formation_volume_factor(P_n, fluid, grid);
bo_n = 1 ./ FVF_n.Bo;
bg_n = 1 ./ FVF_n.Bg;

FVF_n_1 = formation_volume_factor(P_n_1, fluid, grid);
bo_n_1 = 1 ./ FVF_n_1.Bo;
bg_n_1 = 1 ./ FVF_n_1.Bg;

Rs_n = gas_oil_ratio(po_n, fluid, grid);
Rs_n_1 = gas_oil_ratio(po_n_1, fluid, grid);

for i = 1 : grid.blocknums
    %first calculate the accumulation for gas
    accum_vec(2*i-1) = (grid.Vij(i)*rock.porosity/(5.615*delta_t))*((1-sg_n_1(i))...
        *Rs_n_1(i)*bo_n_1(i)+sg_n_1(i)*bg_n_1(i)) - (grid.Vij(i)*rock.porosity/(5.615*delta_t))* ...
        ((1-sg_n(i))*Rs_n(i)*bo_n(i)+sg_n(i)*bg_n(i));
    % then caculate the accumulation for oil
    accum_vec(2*i) = (grid.Vij(i)*rock.porosity/(5.615*delta_t))*((1-sg_n_1(i))*bo_n_1(i)-...
        (1-sg_n(i))*bo_n(i));
end

accumulation = sparse(accum_vec);


end