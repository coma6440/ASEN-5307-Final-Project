%% Housekeeping
close all
clear
clc

%% Load in the data
h_ws = sortrows(readtable("output/HORIZONTAL_WIND_SPEED.csv", 'MissingRule', 'omitrow'));
ws = h_ws.Var5;
t = h_ws.Var1;
h_wd = sortrows(readtable("output/WIND_DIRECTION.csv", 'MissingRule', 'omitrow'));
r_idx = ws < 1;
ws(r_idx) = [];
t(r_idx) = [];
%G = groupsummary(h_ws, 'Var1', 'mean', 'Var5')
%% Fit Weibull Distribution
[p_hat, CI] = wblfit(ws);
pdf_ws = wblpdf(0:0.1:20, p_hat(1), p_hat(2)); 
%% Plotting distribution
figure
histogram(ws, 'Normalization', 'pdf')
hold on
plot(0:0.1:20, pdf_ws, 'LineWidth',2)
xlabel("Wind Speed [m/s]")
ylabel("Probability Density")

%% Plotting variation with moving average
figure
hold on
grid on
scatter(t, ws,'.')
mv_ws =  movmean(ws,200);
plot(t, movmean(ws,200), 'LineWidth', 2)
xlabel("Sol")
ylabel("Wind Speed [m/s]")

%% Time Evolution of Wind Direction
%Histogram not working well due to size of measurements, need to find way
%to remedy
G = groupsummary(h_wd, 'Var1', 'mean', 'Var5');
figure
hold on
scatter(h_wd.Var1, h_wd.Var5,'.')
axis tight
plot(G.Var1, G.mean_Var5)
plot(G.Var1, smooth(G.mean_Var5, 20),'LineWidth', 2)
xlabel("Sol")
ylabel("Wind Direction [Degree]")
legend("Raw Data", "Sol Averaged Data", "Smoothed, Averaged Data")
%Makes sense for wind-drection to be primarily coming from the south,
%curiosity is located in southern hemisphere https://marsed.asu.edu/mep/wind
%Additionally shows significant South EAST !
