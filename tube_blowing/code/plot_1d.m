% Plot 1d integrations

clear; close all; clc;

% folder of datasheets
dir_data = '..\\datasheets\\';
% folder to save in 
dir_save = '..\\1d_select\\';
% file names for relevant data
phi_file = 'phi_1D.csv';
q_file = 'q_1D.csv';
I_phi_file = 'WAXS_IvsPhi.csv';
I_phi_sub_file = 'WAXS_Sub_IvsPhi.csv';
I_q_file = 'WAXS_IvsQ.csv';
I_q_sub_file = 'WAXS_Sub_IvsQ.csv';
% save images?
save = true;
% frames to plot
frames = [0];

% figure parameters
l_fs = 12; % legend font size
phi_min = 0;
phi_max = 180;

%% Load data
phi = csvread(strcat(dir_data, phi_file));
q = csvread(strcat(dir_data, q_file));
I_phi_sub = csvread(strcat(dir_data, I_phi_sub_file));
I_phi = csvread(strcat(dir_data, I_phi_file));
I_q_sub = csvread(strcat(dir_data, I_q_sub_file));
I_q = csvread(strcat(dir_data, I_q_file));

% create labels
labels = {};
for f=frames
    labels = [labels, sprintf('Frame %d', f)];
end

%% Plot IvsQ

% remove NaN
I_q_sub(isnan(I_q_sub)) = 0;

% set axis limits
xmin = 0.68;
xmax = 3.0;
ymin = -100;
ymax = max(max(I_q_sub(:, frames+1)));

figure(3000)
hold on
plot (q, I_q_sub(:, frames+1), 'linewidth', 2)
axis ([xmin xmax ymin ymax])
xlabel ('q(1/Å)')
ylabel ('Intensity[a.u.]')
title('Intensity vs. q (bkgd-subt)')
legend(labels, 'FontSize', l_fs)
set(legend,'location','best')
set(gcf,'PaperPositionMode','auto');
box on
h=findall(gcf,'tag','legend'); % find the handle to the legend
if save
    filename = strcat(dir_save, 'IvsQ1_select');
    print(gcf, filename, '-dpng', '-r300'); 
end

%% Plot I(q) without background subtracted

% remove NaN
I_q(isnan(I_q)) = 0;

% set axis limits
xmin = 0.68;
xmax = 3.0;
ymin = -40;
% 1 includes the background
ymax = max(max(I_q(:, [1, frames+2])));

figure(3500)
hold on
plot (q, horzcat(I_q(:,frames+2), I_q(:,1)), 'linewidth', 2)
axis ([xmin xmax ymin ymax])
xlabel ('q(1/Å)')
ylabel ('Intensity[a.u.]')
title('Raw Data')
legend([labels, {'bkgd'}])
set(legend,'location','best')
set(gcf,'PaperPositionMode','auto');
box on
h=findall(gcf,'tag','legend'); % find the handle to the legend
if save
    filename = strcat(dir_save, 'IvsQ1_raw_select');
    print(gcf, filename, '-dpng', '-r300'); 
end

%% Plot I(phi)

% % remove NaN
I_phi_sub(isnan(I_phi_sub)) = 0;

% set axis limits
xmin = phi_min;
xmax = phi_max;
ymin = -100;
ymax = max(max(I_phi_sub(:, frames+1)));

figure(4000)
plot (phi, I_phi_sub(:, frames+1), 'linewidth', 2)
axis ([xmin xmax ymin ymax])
xlabel ('Phi[°]')
ylabel ('Intensity[a.u.]')
title('Intensity vs. Phi (bkgd-subt)')
legend(labels, 'FontSize', l_fs)
set(legend,'location','best')
set(gcf,'PaperPositionMode','auto');
box on
if save
    filename = strcat(dir_save, 'IvsPhi1_select');
    print(gcf, filename, '-dpng', '-r300'); 
end

%% Plot I(phi) raw (no bkgd subtraction)
% % remove NaN
I_phi(isnan(I_phi)) = 0;

% axis limits
xmin = phi_min;
xmax = phi_max;
ymin = -100;
ymax = max(max(I_phi(:, [1, frames+2])));

figure(4500)
plot (phi, horzcat(I_phi(:,frames+2), I_phi(:,1)), 'linewidth', 2)
axis ([xmin xmax ymin ymax])
xlabel ('Phi[°]')
ylabel ('Intensity[a.u.]')
title('Raw Data')
legend([labels, {'bkgd'}])
set(legend,'location','best')
set(gcf,'PaperPositionMode','auto');
box on
if save
    filename = strcat(dir_save, 'IvsPhi1_raw_select');
    print(gcf, filename, '-dpng', '-r300'); 
end
