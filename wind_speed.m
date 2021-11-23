%% Housekeeping
close all
clear
clc

%% Load in the data
h_ws = sortrows(readtable("output/HORIZONTAL_WIND_SPEED.csv", 'MissingRule', 'omitrow'));
ws = h_ws.Var5;
t = h_ws.Var1;
h_wd = sortrows(readtable("output/WIND_DIRECTION.csv", 'MissingRule', 'omitrow'));
%G = groupsummary(h_ws, 'Var1', 'mean', 'Var5')
%% Fitting a Distribution
pd_wb = fitdist(ws, 'Weibull');
pd_nm = fitdist(ws, 'Normal');
pd_gm = fitdist(ws, 'Gamma');

%% QQ Plot
figure
qqplot(ws)

%% Goodness of fit
[h_nm,p_nm] = chi2gof(ws);
[h_wb,p_wb] = chi2gof(ws,'CDF',pd_wb);
[h_gm,p_gm] = chi2gof(ws,'CDF',pd_gm);
%All fail, perhaps the data sample size is too large?
ws_bs = datasample(ws,1000);
[h_nm_bs,p_nm_bs] = chi2gof(ws_bs);
[h_wb_bs,p_wb_bs] = chi2gof(ws_bs,'CDF',pd_wb);
[h_gm_bs,p_gm_bs] = chi2gof(ws_bs,'CDF',pd_gm);

%% Plotting distribution and visualizing
s_test = 0:0.1:20;
figure
histogram(ws, 'Normalization', 'pdf')
hold on
plot(s_test, pdf(pd_wb,s_test), 'LineWidth',2)
plot(s_test, pdf(pd_nm,s_test), 'LineWidth', 2)
plot(s_test, pdf(pd_gm,s_test), 'LineWidth', 2)
xlabel("Wind Speed [m/s]")
ylabel("Probability Density")
legend("Data", "Weibull Fit", "Normal Fit", "Gamma Fit")


%% Plotting variation with moving average
figure
hold on
grid on
scatter(t, ws,'.')
plot(t, smooth(ws,75), 'LineWidth', 2)
xlabel("Sol")
ylabel("Wind Speed [m/s]")
legend("Data", "Moving average")

%% Time Evolution of Wind Direction
%Histogram not working well due to size of measurements, need to find way
%to remedy. Let's try bootstrapping by sol
G = groupsummary(h_wd, 'Var1', 'mean', 'Var5');
figure
hold on
scatter(h_wd.Var1, h_wd.Var5,'.')
axis tight
plot(G.Var1, G.mean_Var5)
plot(G.Var1, smooth(G.mean_Var5, 75),'LineWidth', 2)
xlabel("Sol")
ylabel("Wind Direction [Degree]")
legend("Raw Data", "Sol Averaged Data", "Smoothed, Averaged Data")
%Makes sense for wind-drection to be primarily coming from the south,
%curiosity is located in southern hemisphere https://marsed.asu.edu/mep/wind
%Additionally shows significant South EAST !