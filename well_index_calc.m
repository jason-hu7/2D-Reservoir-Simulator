function[WI] = well_index_calc(well, rock, grid)
well_loc_x = well.location(1);
well_loc_y = well.location(2);
well_kx = rock.kx(well_loc_x, well_loc_y);
well_ky = rock.ky(well_loc_x, well_loc_y);
d_z = grid.Dbot-grid.Dtop;
% calculate r0
r_0 = 0.28*((well_ky/well_kx)^0.5 * grid.delta_x(well_loc_x, well_loc_y)^2 ... 
    + (well_kx/well_ky)^0.5 * grid.delta_y(well_loc_x, well_loc_y)^2)^0.5...
    / ((well_ky/well_kx)^0.25 + (well_kx/well_ky)^0.25);

% well index
WI = 0.001127*(2*pi*(well_kx*well_ky)^0.5*d_z)/(log(r_0/well.diameter)+well.skin);

end