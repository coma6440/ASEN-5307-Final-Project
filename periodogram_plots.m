function [f, min_peaks, max_peaks] = periodogram_plots(d_min, d_max, n_peaks)
%Compute the periodogram on detrended data
[pxx_min, f_min] = periodogram(detrend(d_min),[],[], 1);
[pxx_max, f_max] = periodogram(detrend(d_max),[],[], 1);

%Find the peak values
[~,m_idx_min] = maxk(pxx_min,n_peaks);
min_peaks = f_min(m_idx_min);

[~,m_idx_max] = maxk(pxx_max,n_peaks);
max_peaks = f_min(m_idx_max);

%Plot the periodograms
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

%Convert the peaks the periods
min_peaks = 1./min_peaks;
max_peaks = 1./max_peaks;
end