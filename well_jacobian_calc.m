function [J_well] = well_jacobian_calc(P_vec, well, fluid, grid)
pressure = P_vec(1:2:size(P_vec));
% fluid = fluid_property_calc(P_vec, fluid, grid);
kro = fluid.k_ro; krg = fluid.k_rg;
bo = fluid.bo; bg = fluid.bg;
mu_o = fluid.oilViscosity; mu_g = fluid.gasViscosity;
Rs = fluid.Rs;

derivatives = derivatives_calc(P_vec, fluid, grid);
dkro_dsg = derivatives.dkro_dsg;
dkrg_dsg = derivatives.dkrg_dsg;
dmu_g_dp = derivatives.dmu_g_dpo;
dbo_dp = derivatives.dbo_dpo;
dbg_dp = derivatives.dbg_dpo;
dRs_dp = derivatives.dRs_dpo;
J_well = zeros(2*grid.blocknums);
% loop through all the wells
for i = 1 : size(well,2)
    well_blocknum = well.blocknum;
    T_wo = well(i).WI*(kro(well_blocknum)*bo(well_blocknum)/mu_o);
    T_wg = well(i).WI*(krg(well_blocknum)*bg(well_blocknum)/mu_g(well_blocknum));
    dT_wo_dp = well(i).WI*kro(well_blocknum)*((1/mu_o)*dbo_dp(well_blocknum));
    dT_wo_dsg = well(i).WI*(bo(well_blocknum)/mu_o)*dkro_dsg(well_blocknum);
    dT_wg_dp = well(i).WI*krg(well_blocknum)*((1/mu_g(well_blocknum))*dbg_dp(well_blocknum) ...
        -(bg(well_blocknum)/(mu_g(well_blocknum)^2))*dmu_g_dp(well_blocknum));
    dT_wg_dsg = well(i).WI*bg(well_blocknum)/mu_g(well_blocknum)*dkrg_dsg(well_blocknum);
    
    if strcmp(well(i).regime, 'BHP')
        dR_o_dp = -T_wo-dT_wo_dp*(pressure(well_blocknum)-well(i).BHP);
        dR_o_dsg = -dT_wo_dsg*(pressure(well_blocknum)-well(i).BHP);
        dR_g_dp = -(T_wg+Rs(well_blocknum)*T_wo)-(dT_wg_dp+Rs(well_blocknum)*dT_wo_dp ...
            +T_wo*dRs_dp(well_blocknum))*(pressure(well_blocknum)-well(i).BHP);
        dR_g_dsg = -(dT_wg_dsg + Rs(well_blocknum)*dT_wo_dsg)*(pressure(well_blocknum)-well(i).BHP);
        % for the gas equations
        J_well(2*well_blocknum-1, 2*well_blocknum-1) = dR_g_dp;
        J_well(2*well_blocknum-1, 2*well_blocknum) = dR_g_dsg;
        % for the oil equations
        J_well(2*well_blocknum, 2*well_blocknum-1) = dR_o_dp;
        J_well(2*well_blocknum, 2*well_blocknum) = dR_o_dsg;
        
    elseif strcmp(well(i).regime, 'RATE')
        dR_o_dp = 0;
        dR_o_dsg = 0;
        dR_g_dp = -well(i).rate*(1/T_wo*dT_wg_dp-T_wg/(T_wo^2)*dT_wo_dp+dRs_dp(well_blocknum));
        dR_g_dsg = -well(i).rate*(1/T_wo*dT_wg_dsg-T_wg/(T_wo^2)*dT_wo_dsg);
        % for the gas equations
        J_well(2*well_blocknum-1, 2*well_blocknum-1) = dR_g_dp;
        J_well(2*well_blocknum-1, 2*well_blocknum) = dR_g_dsg;
        % for the oil equations
        J_well(2*well_blocknum, 2*well_blocknum-1) = dR_o_dp;
        J_well(2*well_blocknum, 2*well_blocknum) = dR_o_dsg;
    else
        fprintf('wrong regime! Please check well specification. \n')
    end
end
J_well = sparse(J_well);
end