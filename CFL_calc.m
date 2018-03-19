function [simulation] = CFL_calc(well,rock,simulation,grid)
x = well(1).location(1); y = well(1).location(2);
q_t = well(1).gas_rate/178.10760067 + well(1).oil_rate;
PV = grid.Vij(y,x)*rock.porosity/5.615;
simulation.CFL = q_t * simulation.time_step / PV;
end