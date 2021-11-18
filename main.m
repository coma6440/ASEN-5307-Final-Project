%% Housekeeping
close all
clear
clc
%%
% A = readtable('output/WIND_DIRECTION.csv');
% B = readtable('output/HORIZONTAL_WIND_SPEED.csv');
% % figure
% % scatter(A.Var1, A.Var5)
% % polarhistogram(A.Var5, 'DisplayStyle', 'stairs')
% 
% histogram(B.Var5)
% scatter(B.Var2, B.Var5)

%% Load in the Data
temperature = sortrows(readtable("output/temp_stats.csv", 'MissingRule', 'omitrow'));
pressure = sortrows(readtable("output/pressure_stats.csv", 'MissingRule', 'omitrow'));
humidity = sortrows(readtable("output/humidity_stats.csv", 'MissingRule', 'omitrow'));

%% Initial Plots of Results
% plot_summary(temperature, "Daily Temperature Summary");
% plot_summary(pressure, "Daily Pressure Summary");
% plot_summary(humidity, "Daily Humidity Summary");

%% Remove Outliers
temperature = rmoutliers(temperature);
pressure = rmoutliers(pressure);
humidity = rmoutliers(humidity);
%% Interpolation of missing values
% temperature = interp_table(temperature,  temperature.SOL, idx);
% pressure = interp_table(pressure,  pressure.SOL, idx);
% humidity = interp_table(humidity,  humidity.SOL, idx);

%% Fitting Fourier Series
idx = 1:max(temperature.SOL);

% t_s = temperature.SOL;
% t_min = temperature.MINIMUM;
% t_max = temperature.MAXIMUM;
% 
% p_s = pressure.SOL;
% p_min = pressure.MINIMUM;
% p_max = pressure.MAXIMUM;
% 
% h_s = humidity.SOL;
% h_min = humidity.MINIMUM;
% h_max = humidity.MAXIMUM;

% Made using cftool
load('fits.mat');

[t_min, t_max, ~] = plot_fit(idx, temperature, t_min_fit, t_max_fit, "Temperature Fitted Model");
[p_min, p_max, ~] = plot_fit(idx, pressure, p_min_fit, p_max_fit, "Pressure Fitted Model");
[h_min, h_max, ~] = plot_fit(idx, humidity, h_min_fit, h_max_fit, "Humidity Fitted Model");

resid_histogram(t_min_out, t_max_out, "Tempearature Fit Residuals");
resid_histogram(p_min_out, p_max_out, "Pressure Fit Residuals");
resid_histogram(h_min_out, h_max_out, "Humidity Fit Residuals");

%Making table of GOF statistics
gof(1) = t_min_good;
gof(2) = t_max_good;
gof(3) = p_min_good;
gof(4) = p_max_good;
gof(5) = h_min_good;
gof(6) = h_max_good;

%% Continuous Wavelet Transform
%Goal here is to find temporal signals


%% Plots After Further Cleaning
% plot_summary(temperature, "Daily Temperature Summary");
% plot_summary(pressure, "Daily Pressure Summary");
% plot_summary(humidity, "Daily Humidity Summary");
%% Helper Functions
function f = plot_summary(A, title_str)
f = figure;
nexttile
scatter(A.SOL, A.MINIMUM, '.')
xlabel("SOL")
ylabel("MINIMUM")
nexttile
scatter(A.SOL, A.MAXIMUM, '.')
xlabel("SOL")
ylabel("MAXIMUM")
nexttile
scatter(A.SOL, A.MEAN, '.')
xlabel("SOL")
ylabel("MEAN")
nexttile
scatter(A.SOL, A.MEDIAN, '.')
xlabel("SOL")
ylabel("MEDIAN")
nexttile
scatter(A.SOL, A.MODE, '.')
xlabel("SOL")
ylabel("MODE")
nexttile
scatter(A.SOL, A.STD, '.')
xlabel("SOL")
ylabel("STANDARD DEVIATION")
sgtitle(title_str)
end

function [min_fit, max_fit, f] = plot_fit(t, data, min_fit, max_fit, title_str)
min_fit = feval(min_fit,t);
max_fit = feval(max_fit,t);
f = figure;
nexttile
scatter(data.SOL, data.MINIMUM, '.');
hold on
plot(t, min_fit, 'LineWidth', 2)
xlabel("SOL")
ylabel("VALUE")
nexttile
scatter(data.SOL, data.MAXIMUM, '.');
hold on
plot(t, max_fit, 'LineWidth', 2)
xlabel("SOL")
ylabel("VALUE")
legend("Data", "Fitted Model")
sgtitle(title_str)
end

function f = resid_histogram(min_out, max_out, title_str)
f = figure;
nexttile
histogram(min_out.residuals, 'Normalization', 'probability')
nexttile
histogram(max_out.residuals, 'Normalization', 'probability')
sgtitle(title_str)
end