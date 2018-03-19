%% function that performs automatic time stepping
function[delta_t_1] = time_stepping(P_n_1, P_n, simulation)
po_n_k = P_n_1(1:2:size(P_n_1));
po_n = P_n(1:2:size(P_n));
sat_n_k =P_n_1(2:2:size(P_n_1));
sat_n = P_n(2:2:size(P_n));
delta_t_n = simulation.time_step;

omega = 0.1;    % tuning parameter
eta_p = 50;   % maximum desired change for pressure
eta_s = 0.05;   % maximum desired change for gas saturation

min_po = min((1+omega)*eta_p./(abs(po_n_k-po_n)+omega*eta_p));
min_sat = min((1+omega)*eta_s./(abs(sat_n_k-sat_n)+omega*eta_s));
min_x = min(min_po, min_sat);
delta_t_1 = delta_t_n * min_x;

if delta_t_1 > simulation.max_step
    delta_t_1 = simulation.max_step;
end

if (delta_t_1 + simulation.cur_time) > simulation.tot_time
    delta_t_1 = simulation.tot_time - simulation.cur_time;
end

end