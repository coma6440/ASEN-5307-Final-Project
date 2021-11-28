%% Housekeeping
close all
clear
clc

%% Load in the Data
temperature = sortrows(readtable("output/temp_stats.csv", 'MissingRule', 'omitrow'));
pressure = sortrows(readtable("output/pressure_stats.csv", 'MissingRule', 'omitrow'));
humidity = sortrows(readtable("output/humidity_stats.csv", 'MissingRule', 'omitrow'));

%% Initial Plots of Results
f = plot_summary(temperature);
saveas(f, "images/t_summary_init.png")
f = plot_summary(pressure);
saveas(f, "images/p_summary_init.png")
f = plot_summary(humidity);
saveas(f, "images/h_summary_init.png")
%% Remove Outliers
temperature = rmoutliers(temperature);
pressure = rmoutliers(pressure);
humidity = rmoutliers(humidity);

%% Fitting Fourier Series
idx = 1:max(temperature.SOL);

t_min = temperature.MINIMUM;
t_max = temperature.MAXIMUM;
 
p_min = pressure.MINIMUM;
p_max = pressure.MAXIMUM;

h_min = humidity.MINIMUM;
h_max = humidity.MAXIMUM;

% Made using cftool
load('fits.mat');

[~, ~, f] = plot_fit(idx, temperature, t_min_fit, t_max_fit, "Temperature [K]");
saveas(f, "images/t_fit.png")
[~, ~, f] = plot_fit(idx, pressure, p_min_fit, p_max_fit, "Pressure [Pa]");
saveas(f, "images/p_fit.png")
[~, ~, f] = plot_fit(idx, humidity, h_min_fit, h_max_fit, "Humidity [%]");
saveas(f, "images/h_fit.png")

%Making table of GOF statistics
gof(1) = t_min_good;
gof(2) = t_max_good;
gof(3) = p_min_good;
gof(4) = p_max_good;
gof(5) = h_min_good;
gof(6) = h_max_good;

%% Continuous Wavelet Transform
%Goal here is to find temporal signals
f = cwt_plots(idx, t_min, t_max);
saveas(f, "images/t_cwt.png")
f = cwt_plots(idx, p_min, p_max);
saveas(f, "images/p_cwt.png")
f = cwt_plots(idx, h_min, h_max);
saveas(f, "images/h_cwt.png")

%% Periodogram to find peaks
[f, t_peaks_min, t_peaks_max] = periodogram_plots(t_min, t_max, 1);
saveas(f, "images/t_periodogram.png")
[f, p_peaks_min, p_peaks_max] = periodogram_plots(p_min, p_min, 2);
saveas(f, "images/p_periodogram.png")
[f, h_peaks_min, h_peaks_max] = periodogram_plots(h_min, h_max, 1);
saveas(f, "images/h_periodogram.png")

%% Correlation between variables
idx = intersect(intersect(temperature.SOL, pressure.SOL), humidity.SOL);
t_idx = zeros(length(idx),1);
p_idx = zeros(length(idx),1);
h_idx = zeros(length(idx),1);

for i = 1:length(idx)
    t_idx(i) = find(temperature.SOL == idx(i));
    p_idx(i) = find(pressure.SOL == idx(i));
    h_idx(i) = find(humidity.SOL == idx(i));
end

data = [t_min(t_idx), t_max(t_idx), p_min(p_idx), p_max(p_idx), h_min(h_idx), h_max(h_idx)];
R = corrcoef(data);
tick_labels = ["Min(T)","Max(T)","Min(P)","Max(P)","Min(H)","Max(H)"];
f = figure('Visible','off');
h = heatmap(R);
ax = gca;
ax.XDisplayLabels = tick_labels;
ax.YDisplayLabels = tick_labels;
saveas(f, "images/var_corr.png")