%% calculate the fluid properties based on the pressure and saturation given
function [fluid] = fluid_property_calc(P_vec, fluid, grid)

% gas viscosity is a function of pressure
gasViscosity = @(pressure) 3E-10 * pressure.^2 + 1E-6 * pressure + 0.0133;

% Relative Permeability
k_ro = @(s_g) (1 - s_g) .^ 1.5;
k_rg = @(s_g) s_g .^ 2.0;

pressure = P_vec(1:2:size(P_vec));
sg = P_vec(2:2:size(P_vec));

FVF = formation_volume_factor(P_vec, fluid, grid);
fluid.bo = 1 ./ FVF.Bo;
fluid.bg = 1 ./ FVF.Bg;
fluid.gasViscosity  = gasViscosity(pressure);
fluid.k_ro = k_ro(sg);
fluid.k_rg = k_rg(sg);
fluid.Rs = gas_oil_ratio(pressure, fluid, grid);

end