%% provide input for the simulation
function [grid, rock, fluid, well, simulation] = simulation_input_4_2()

% Reservoir grid description ----------------------------------------------
grid.Nx = 23;          % number of grid blocks in x direction
grid.Ny = 23;          % number of grid blocks in y direction
grid.Nz = 1;           % number of grid blocks in z direction
grid.blocknums = grid.Nx * grid.Ny * grid.Nz;

grid.Lx = 6900;        % x length (ft) of the reservoir
grid.Ly = 6900;        % y length (ft) of the reservoir
grid.Lz = 100;         % z length (ft) of the reservoir
grid.Dtop = 1200;      % top depth (ft) of the reservoir
grid.Dbot = 1300;    % bottom depth (ft) of the reservoir

delta_x = repmat(repmat(grid.Lx / grid.Nx, 1, grid.Nx), grid.Ny, 1);
delta_y = repmat(repmat(grid.Ly / grid.Ny, 1, grid.Nx), grid.Ny, 1);
delta_z = grid.Lz / grid.Nz;
grid.delta_x = repmat(grid.Lx / grid.Nx, grid.Ny, grid.Nx);
grid.delta_y = repmat(grid.Ly / grid.Ny, grid.Ny, grid.Nx);
grid.Vij = delta_x .* delta_y .* delta_z;

% grid.wellLoc_1 = [12,12];% location of the well in the base case
% grid.wellLoc_2 = [5,5];% location of the well in the base case
% grid.inspectLoc_1 = [18,18];
% grid.inspectLoc_2 = [5,5];

% Rock properties ---------------------------------------------------------
kx_mat = 10* ones(grid.Ny,grid.Nx);  
kx_mat(12,:) = 200;
ky_mat = 10* ones(grid.Ny,grid.Nx);  
ky_mat(:,12) = 200;
rock.kx = kx_mat;         % permeability (md) in x direction
rock.ky = ky_mat;         % permeability (md) in y direction

rock.porosity = 0.22;  % porosity
rock.CR = 0;           % rock compressibility

% Fluid properties --------------------------------------------------------
fluid.oilCf = 0.8E-5;       % oil compressibility (psi^-1)
fluid.oilDensity = 49.1;    % oil density at surface conditions(lbm/ft^3)
fluid.gasDensity = 0.06055; % gas density at surface conditions(lbm/scf)
fluid.oilViscosity = 2.5;   % oil viscosity (cp)
fluid.bubblePressure = 3400;    % bubble point pressure (psi)
fluid.surfPressure = 14.7;   % surface pressure (psi)
fluid.pinit = 4500;         % intial pressure of the reservoir at 4500 psi @ 1250ft

% Well properties----------------------------------------------------------
well(1).diameter = 0.5; % well diameter is 0.5ft
well(1).location = [12,12];    % well location
well(1).blocknum = (well(1).location(2)-1)*grid.Ny + well(1).location(1);
% well(1).regime = 'BHP'; % This is a fixed BHP well
% well(1).BHP = 6000;     % BHP fixed regime at 2000 psi
well(1).regime = 'RATE';
well(1).rate = -400;
well(1).min_BHP = 2000;
well(1).skin = 0;
% calculate well index
well(1).WI = well_index_calc(well(1), rock, grid);


% well(2).diameter = 0.5; % well diameter is 0.5ft
% well(2).location = [8,8];    % well location
% well(2).BHP = 2000;     % BHP fixed regime at 2000 psi

% simulation parameters ---------------------------------------------------
simulation.time_step = 0.01;    % (days) of simulation
simulation.max_step = 20;  % max time step size
simulation.tot_time = 1500;
simulation.cur_time = 0;



end