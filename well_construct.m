function[W_vec, well] = well_construct(P_n, well, grid, fluid)

pressure = P_n(1:2:size(P_n));
R_o = zeros(2*grid.blocknums, 1);
R_g = zeros(2*grid.blocknums, 1);

i = 1;
% loop through all the wells
while i <= size(well,2)
    well_blocknum = well(i).blocknum;
    T_wo = well(i).WI*(fluid.k_ro(well_blocknum)*fluid.bo(well_blocknum)/fluid.oilViscosity);
    T_wg = well(i).WI*(fluid.k_rg(well_blocknum)*fluid.bg(well_blocknum)/fluid.gasViscosity(well_blocknum));
    
    if strcmp(well(i).regime, 'BHP')
        R_g(well_blocknum*2-1) = -(T_wg+fluid.Rs(well_blocknum)*T_wo)*(pressure(well_blocknum)-well(i).BHP);
        R_o(well_blocknum*2) = -T_wo * (pressure(well_blocknum) - well(i).BHP); 
        well(i).oil_rate = -R_o(well_blocknum*2);
        well(i).gas_rate = -R_g(well_blocknum*2-1) / 1000 * 5.615;
        i = i + 1;
    elseif strcmp(well(i).regime, 'RATE')
        R_g(well_blocknum*2-1) = -(T_wg/T_wo+fluid.Rs(well_blocknum))*well(i).rate;
        R_o(well_blocknum*2) = -well(i).rate;
        well(i).BHP = (-well(i).rate+T_wo*pressure(well_blocknum))/T_wo;
        if well(i).BHP < well(i).min_BHP
            well(i).regime = 'BHP';
            well(i).BHP = well(i).min_BHP;
            fprintf('switch rate control regime to pressure control regime and repeat last step\n');
        else
            well(i).oil_rate = -R_o(well_blocknum*2);
            well(i).gas_rate = -R_g(well_blocknum*2-1) / 1000 * 5.615;
            i = i + 1;
        end
    else
        fprintf('wrong regime! Please check well specification. \n')
    end
%     well(i).oil_rate = -R_o(well_blocknum*2);
%     well(i).gas_rate = -R_g(well_blocknum*2-1) / 1000 * 5.615;
end
W_vec = R_o + R_g;
W_vec = sparse(W_vec);


end