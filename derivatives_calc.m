%% Calculate the derviatives of bo, bg, Rs, kro and mu_g at the center of each blocks
function [derivatives] = derivatives_calc(P_vec, fluid, grid)

p_atm = fluid.surfPressure;
p = P_vec(1:2:size(P_vec));
sg = P_vec(2:2:size(P_vec));

dbo_dpo = zeros(grid.blocknums, 1);
dbg_dpo = zeros(grid.blocknums, 1);
dRs_dpo = zeros(grid.blocknums, 1);

p_bub = fluid.bubblePressure;
co = fluid.oilCf;


%% calculate the derivatives of kro with respect to sg
dkro_dsg = -1.5.*((1-sg).^0.5);

%% calculate the derivatives of krg with respect to sg
dkrg_dsg = 2 .* sg;

%% calculate the derivatives of mu_g with respect to sg
dmu_g_dpo = 6E-10 .* p + 1E-6;

%% calculate the derivatives of bo with respect to po
for i = 1 : grid.blocknums
    if p(i) < p_bub
        dbo_dpo(i) = -8.0E-5 * exp(8.0E-5 * (p_atm - p(i)));
        dbg_dpo(i) = 1.7E-3 * exp(1.7E-3 * p(i) - 1.7E-3 * p_atm);
        dRs_dpo(i) = (178.11^2 / 5.615) * 1.3 * (p(i)^0.3) / (p_bub ^ 1.3);
    else
        dbo_dpo(i) = exp(8E-5*(p_atm-p_bub))*exp(co*(p(i)-p_bub))*co;
        dbg_dpo(i) = 0;
        dRs_dpo(i) = 0;
    end
end

derivatives.dkro_dsg = dkro_dsg;
derivatives.dkrg_dsg = dkrg_dsg;
derivatives.dmu_g_dpo = dmu_g_dpo;
derivatives.dbo_dpo = dbo_dpo;
derivatives.dbg_dpo = dbg_dpo;
derivatives.dRs_dpo = dRs_dpo;

end