%% check on the convergence criteria for the Newton loop
function [check] = convergence_check(Residual, P_n_k, P_n, grid, rock, fluid, simulation)

delta_t = simulation.time_step;
norm_vec = zeros(2*grid.blocknums,1);
po_n_k = P_n_k(1:2:size(P_n_k));
po_n = P_n(1:2:size(P_n));
sat_n_k =P_n_k(2:2:size(P_n_k));
sat_n = P_n(2:2:size(P_n));
check = 0;

for i = 1 : grid.blocknums
    PV_i = rock.porosity * grid.Vij(i);
    % inf norm for gas
    norm_vec(2*i-1) = abs(5.615*delta_t*Residual(2*i-1)/fluid.bg(i)/PV_i);
    % inf norm for oil
    norm_vec(2*i) = abs(5.615*delta_t*Residual(2*i)/fluid.bo(i)/PV_i);
end
inf_norm = max(norm_vec);
max_sat_change = max(abs(sat_n_k - sat_n));
max_po_change = max(abs((po_n_k - po_n) / mean(po_n_k)));

e_1 = 10^-3; e_2 = 0.01; e_3 = 0.001;

if (inf_norm <= e_1) & (max_sat_change <= e_2) & (max_po_change <= e_3)
    check = 1;
end


end