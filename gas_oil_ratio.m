%% calculate the gas oil ratio given oil pressure
function [Rs] = gas_oil_ratio(pressure, fluid, grid)

Rs = zeros(grid.blocknums,1);
for i = 1 : grid.blocknums
    if pressure(i) <= fluid.bubblePressure
        Rs(i) = (178.11^2)/5.615 * ((pressure(i)/fluid.bubblePressure)^1.3);
    else
        Rs(i) = (178.11^2)/5.615;
    end  
end

end
