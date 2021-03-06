%% Housekeeping
close all
clear
clc

%% Load in the data
h_ws = sortrows(readtable("output/HORIZONTAL_WIND_SPEED.csv", 'MissingRule', 'omitrow'));
ws = h_ws.Var5;
t = h_ws.Var1;
h_wd = sortrows(readtable("output/WIND_DIRECTION.csv", 'MissingRule', 'omitrow'));

%% Fitting Theoretical Distribution
pd_wb = fitdist(ws, 'Weibull');
pd_nm = fitdist(ws, 'Normal');
pd_gm = fitdist(ws, 'Gamma');

%% QQ Plot
f = figure('Visible','off');
qqplot(ws)
saveas(f, "images/ws_qq.png")

%% Goodness of fit
[h_nm,p_nm] = chi2gof(ws);
[h_wb,p_wb] = chi2gof(ws,'CDF',pd_wb);
[h_gm,p_gm] = chi2gof(ws,'CDF',pd_gm);

%Resample data with smaller sample size
ws_bs = datasample(ws,1000);

%Try tests again
[h_nm_bs,p_nm_bs] = chi2gof(ws_bs);
[h_wb_bs,p_wb_bs] = chi2gof(ws_bs,'CDF',pd_wb);
[h_gm_bs,p_gm_bs] = chi2gof(ws_bs,'CDF',pd_gm);

%% Plotting distribution and visualizing theoretical dists
%Wind speeds for testing
s_test = 0:0.1:20;
f = figure('Visible','off');
histogram(ws, 'Normalization', 'pdf')
hold on
plot(s_test, pdf(pd_wb,s_test), 'LineWidth',2)
plot(s_test, pdf(pd_nm,s_test), 'LineWidth', 2)
plot(s_test, pdf(pd_gm,s_test), 'LineWidth', 2)
xlabel("Wind Speed [m/s]")
ylabel("Probability Density")
legend("Data", "Weibull Fit", "Normal Fit", "Gamma Fit")
saveas(f, "images/ws_fitted.png")

%% Plotting variation with moving average
f = figure('Visible','off');
hold on
grid on
scatter(t, ws,'.')
plot(t, smooth(ws,100), 'LineWidth', 2)
xlabel("Sol")
ylabel("Wind Speed [m/s]")
legend("Data", "Moving average")
saveas(f, "images/ws_time.png")

%% Get the daily time of the wind directions
sols = h_wd.Var1;
t_wd = h_wd.Var2;
wd = h_wd.Var5;
n = length(h_wd.Var1);

if isfile('wd_time.mat')
    load('wd_time.mat');
else
    for i = 1:n
        t_str = h_wd.Var3{i};
        time(i) = duration(t_str(7:end),'Format','hh:mm:ss');
    end
    save('wd_time.mat', time);
end

%% Find the angular mean of wind direction for each day
G = groupsummary(h_wd, 'Var1', @(x) ang_mean(x), 'Var5');
t_mean = G.Var1;
wd_mean = rad2deg(mod(G.fun1_Var5,2*pi));

%% Plot the time series for wind direction
f = figure('Visible','off');
scatter(sols,wd,'.')
axis tight
xlabel("Sol")
ylabel("Wind Direction [Degrees]")
hold on
scatter(t_mean, wd_mean,'.')
legend(["Data", "Daily Angular Mean"])
saveas(f, "images/wd_time.png")

%% Histogram of observation times

f = figure('Visible','off');
histogram(time,20)
xlabel("LMST")
ylabel("Observations")
axis tight
saveas(f, "images/wd_obs_times.png")
%Majority of samples taken during the day

%% Polar histogram of observations
f = figure('Visible','off');
polarhistogram(deg2rad(wd))
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
angles = 0:45:360;
pax.ThetaTick = angles;
labels = {'N','NE','E','SE','S','SW','W','NW'};
pax.ThetaTickLabel = labels;
saveas(f, "images/wd_histogram.png")

%% Splitting Wind-Direction by Time
t_rise = duration(hours(6));
t_fall = duration(hours(18));
day = time > t_rise & time < t_fall;
night = time < t_rise | time > t_fall;

% Daytime
f = figure('Visible','off');
polarhistogram(deg2rad(wd(day)))
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
angles = 0:45:360;
pax.ThetaTick = angles;
labels = {'N','NE','E','SE','S','SW','W','NW'};
pax.ThetaTickLabel = labels;
saveas(f, "images/wd_histogram_day.png")

%Nighttime
f = figure('Visible','off');
polarhistogram(deg2rad(wd(night)))
pax = gca;
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
angles = 0:45:360;
pax.ThetaTick = angles;
labels = {'N','NE','E','SE','S','SW','W','NW'};
pax.ThetaTickLabel = labels;
saveas(f, "images/wd_histogram_night.png")