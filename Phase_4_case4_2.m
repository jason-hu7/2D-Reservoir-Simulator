clear; clc
tic
% Phase 3 of the ENERGY 223 project
% In this phase of the project we wil build a simulator for single-phase (oil)
% solution will be implicit but evaluate the transmissibilities at the previous time. 

%% Input parameters
[grid, rock, fluid, well, simulation] = simulation_input_4_2();

%% Initialzation
% initial primary variable vectors
rock_geo_trans = rock_geo_transmissbility(grid, rock);
% initial p vector and the first guess for newton's iterations
P_n = zeros(2*grid.blocknums,1);
for i = 1 : 2*grid.blocknums
    if mod(i,2) == 0
        P_n(i) = 0.4;
    else
        P_n(i) = fluid.pinit;
    end
end

% all the simulation output we want to observe
init_pressure = P_n(1:2:size(P_n));
BHP_hist = fluid.pinit;
BPR_hist = init_pressure(well(1).blocknum);
FPR_hist = mean(init_pressure);
time_hist = 0;
FOPR_hist = 0;
FGPR_hist = 0;
CFL_hist = 0;
newton_hist = 0;


%% main body of the simulation
t = simulation.cur_time;
iter = 0;
while t < simulation.tot_time
    iter = iter + 1;
    fprintf('current time is %d days \n',t);
    [P_n_1,simulation, well] = nonlinear_solver(P_n, fluid, rock, grid, rock_geo_trans, simulation, well);
    
    % record BHP, well block pressure and average reservoir pressure
    pressure_n_1 = P_n_1(1:2:size(P_n_1));
    sat_n_1 = P_n_1(2:2:size(P_n_1));
    
    FPR_hist(iter+1) = mean(pressure_n_1);
    BPR_hist(iter+1) = pressure_n_1(well(1).blocknum);
    BHP_hist(iter+1) = well(1).BHP;
    FOPR_hist(iter+1) = -well(1).oil_rate;
    FGPR_hist(iter+1) = -well(1).gas_rate;
    simulation = CFL_calc(well,rock,simulation,grid);
    CFL_hist(iter+1) = -simulation.CFL;
    newton_hist(iter+1) = simulation.newton_num;

    t = t + simulation.time_step;
    simulation.cur_time = t;
    time_hist(iter+1) = t;
    
    if t > 100 & t <300
        well(1).rate = 0;
    elseif t > 500 & t < 700
        well(1).rate = 0;
    elseif t > 900 & t < 1100
        well(1).rate = 0;
    elseif t > 1300 & t < 1500
        well(1).rate = 0;
    else
        well(1).rate = -400;
    end
    
    delta_t_1 = time_stepping(P_n_1, P_n, simulation);
    simulation.time_step = delta_t_1;
    fprintf('The current time step is %d \n', simulation.time_step);
    
    P_n = P_n_1;
    
end

%% load eclipse data
cur_path = pwd(); eclipse_path = fullfile(cur_path, 'Eclipse', 'PHASE4_1.RSM');
eclipse_data = importdata(eclipse_path, ' ', 10);

figure(1)
pressure_field = reshape(P_n_1(1:2:size(P_n_1)), 23, 23)';
surf(pressure_field);
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
xlabel('reservoir length (m)'); xticklabels = 0:1500:6900;
xticks = linspace(1, size(pressure_field, 2), numel(xticklabels));
ylabel('reservoir length (m)'); yticklabels = 0:1500:6900;
yticks = linspace(1, size(pressure_field, 1), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
c = colorbar; set(c,'FontSize',14);
colormap jet; colorTitleHandle = get(colorbar,'Title');
titleString = 'pressure (psi)'; set(colorTitleHandle ,'String',titleString);
saveas(gcf, 'pressure field case 4_2.png')

figure(2)
sat_field = reshape(P_n_1(2:2:size(P_n_1)), 23, 23)';
surf(sat_field);
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
xlabel('reservoir length (m)'); xticklabels = 0:1500:6900;
xticks = linspace(1, size(pressure_field, 2), numel(xticklabels));
ylabel('reservoir length (m)'); yticklabels = 0:1500:6900;
yticks = linspace(1, size(pressure_field, 1), numel(yticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels); 
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);
c = colorbar; set(c,'FontSize',14);
colormap jet; colorTitleHandle = get(colorbar,'Title');
titleString = 'saturation (%)'; set(colorTitleHandle ,'String',titleString);
saveas(gcf, 'sat field case 4_2.png')

figure(3)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,4), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, FOPR_hist, 'DisplayName', 'Jason Simulator');
xlabel('time (days)');ylabel('oil production rate(bbl/day)');
title('Oil Rate History'); legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'oil rate history case 4_2.png')

figure(4)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,3), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, FGPR_hist, 'DisplayName', 'Jason Simulator'); 
xlabel('time (days)');ylabel('Gas production rate(MSCF/day)');
title('Gas Rate Hisotry');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'gas rate history case 4_2.png')

figure(5)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,8), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, FPR_hist, 'DisplayName', 'Jason Simulator'); 
xlabel('time (days)');ylabel('average reservoir pressure(psi)');
title('Average Reservoir Pressure');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'FPR case 4_2.png')

figure(6)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,10), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, BPR_hist, 'DisplayName', 'Jason Simulator');
xlabel('time (days)');ylabel('well block pressure(psi)');
title('Well Block Pressure');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'BPR case 4_2.png')

figure(7)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,9), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, BHP_hist, 'DisplayName', 'Jason Simulator'); 
xlabel('time (days)');ylabel('BHP(psi)');
title('Well Bottom Hole Pressure');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'BHP case 4_2.png')

figure(8)
% plot(eclipse_data.data(:,1), eclipse_data.data(:,9), 'DisplayName', 'Eclipse'); hold on;
plot(time_hist, newton_hist, 'DisplayName', 'Nx=Ny=23'); hold on;
title('Newton loop number at each time step');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'Newton steps case 4_2.png')

figure(9)
plot(time_hist, CFL_hist, 'DisplayName', 'Jason Simulator');
title('max CFL number');legend('show')
set(gcf, 'Color', [1,1,1]); % white background
set(gca,'FontSize',14); % bigger font size
saveas(gcf, 'CFL case 4_2.png')


toc

