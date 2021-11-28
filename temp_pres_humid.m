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

f = resid_histogram(t_min_out, t_max_out);
saveas(f, "images/t_resid_hist.png")
f = resid_histogram(p_min_out, p_max_out);
saveas(f, "images/p_resid_hist.png")
f = resid_histogram(h_min_out, h_max_out);
saveas(f, "images/h_resid_hist.png")


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
%% Helper Functions
function f = plot_summary(A)
A = rmoutliers(A);
f = figure('Visible','off');
nexttile
scatter(A.SOL, A.MINIMUM, '.')
xlabel("Sol")
ylabel("MINIMUM")
nexttile
scatter(A.SOL, A.MAXIMUM, '.')
xlabel("Sol")
ylabel("MAXIMUM")
nexttile
scatter(A.SOL, A.MEAN, '.')
xlabel("Sol")
ylabel("MEAN")
nexttile
scatter(A.SOL, A.MEDIAN, '.')
xlabel("Sol")
ylabel("MEDIAN")
nexttile
scatter(A.SOL, A.MODE, '.')
xlabel("Sol")
ylabel("MODE")
nexttile
scatter(A.SOL, A.STD, '.')
xlabel("Sol")
ylabel("STANDARD DEVIATION")
end

function [min_fit, max_fit, f] = plot_fit(t, data, min_fit, max_fit, yaxis_str)
min_fit = feval(min_fit,t);
max_fit = feval(max_fit,t);
f = figure('Visible','off');
nexttile
scatter(data.SOL, data.MINIMUM, '.');
hold on
plot(t, min_fit, 'LineWidth', 2)
xlabel("Sol")
ylabel(yaxis_str)
axis tight
title("Minimum")

nexttile
scatter(data.SOL, data.MAXIMUM, '.');
hold on
plot(t, max_fit, 'LineWidth', 2)
xlabel("Sol")
ylabel(yaxis_str)
legend("Data", "Fitted Model")
axis tight
title("Maximum")
end

function f = resid_histogram(min_out, max_out)
f = figure('Visible','off');
nexttile
histogram(min_out.residuals)
ylabel("Counts")
xlabel("Residual Value")
title("Minimum")
nexttile
histogram(max_out.residuals)
ylabel("Counts")
xlabel("Residual Value")
title("Maximum")
end

function f = cwt_plots(t, d_min, d_max)
[cmin,fmin] = cwt(detrend(d_min));
[cmax,fmax] = cwt(detrend(d_max));

f = figure('Visible','off');
tiledlayout(1,2)

nexttile
image('XData',t,'YData',fmin,'CData', abs(cmin),'CDataMapping','scaled')
set(gca,'YScale', 'log')
xlabel("Sol")
ylabel("Cycles/Sol")
axis tight
title("Minimum")

nexttile
image('XData',t,'YData',fmax,'CData', abs(cmax),'CDataMapping','scaled')
set(gca,'YScale', 'log')
xlabel("Sol")
ylabel("Cycles/Sol")
axis tight
title('Maximum')

end

function [f, min_peaks, max_peaks] = periodogram_plots(d_min, d_max, n_peaks)
[pxx_min, f_min] = periodogram(detrend(d_min),[],[], 1);
[pxx_max, f_max] = periodogram(detrend(d_max),[],[], 1);

[~,m_idx_min] = maxk(pxx_min,n_peaks);
min_peaks = f_min(m_idx_min);

[~,m_idx_max] = maxk(pxx_max,n_peaks);
max_peaks = f_min(m_idx_max);

f = figure('Visible','off');
tiledlayout(1,2)
nexttile
semilogx(f_min, 20*log10(pxx_min));
hold on
scatter(min_peaks, 20*log10(pxx_min(m_idx_min)))
ylabel("Power [dB/(Cycles/Sol)]")
xlabel("Frequency [Cycles/Sol]")
title("Minimum")
grid on

nexttile
semilogx(f_max, 20*log10(pxx_max));
hold on
scatter(max_peaks, 20*log10(pxx_max(m_idx_max)))
ylabel("Power [dB/(Cycles/Sol)]")
xlabel("Frequency [Cycles/Sol]")
title("Maximum")
grid on
min_peaks = 1./min_peaks;
max_peaks = 1./max_peaks;
end