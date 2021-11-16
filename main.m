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
idx = 1:max(temperature.SOL);

% temperature = interp_table(temperature,  temperature.SOL, idx);
% pressure = interp_table(pressure,  pressure.SOL, idx);
% humidity = interp_table(humidity,  humidity.SOL, idx);

% Can try filtering or least squares fitting?

%% Plots After Further Cleaning
plot_summary(temperature, "Daily Temperature Summary");
plot_summary(pressure, "Daily Pressure Summary");
plot_summary(humidity, "Daily Humidity Summary");
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
function out = interp_table(tab, a_idx, f_idx)
SOL = f_idx';
MINIMUM = interp1(a_idx, tab.MINIMUM, f_idx, 'linear')';
MAXIMUM = interp1(a_idx, tab.MAXIMUM, f_idx, 'linear')';
MEDIAN = interp1(a_idx, tab.MEDIAN, f_idx, 'linear')';
MEAN = interp1(a_idx, tab.MEAN, f_idx, 'linear')';
MODE = interp1(a_idx, tab.MODE, f_idx, 'linear')';
STD = interp1(a_idx, tab.STD, f_idx, 'linear')';
out = table(SOL, MAXIMUM, MINIMUM, MEDIAN, MEAN, MODE, STD);
end