function [FVF] = formation_volume_factor(P, fluid, grid)

pressure = P(1:2:size(P));
Dp = zeros(grid.blocknums,1);
FVF.Bo = zeros(grid.blocknums,1);

for i = 1 : grid.blocknums
    if pressure(i) < fluid.bubblePressure
        Dp(i) = fluid.surfPressure - pressure(i);
        FVF.Bo(i) = exp(-8.0E-5 * Dp(i));
    else
        Dp(i) = fluid.surfPressure - fluid.bubblePressure;
        FVF.Bo(i) = exp(-8.0E-5 * Dp(i)) * exp(-fluid.oilCf * (pressure(i) - fluid.bubblePressure));
    end
end

FVF.Bg = exp(1.7E-3 * Dp);

end