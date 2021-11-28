function f = cwt_plots(t, d_min, d_max)
%Compute the CWT on detrended data
[cmin,fmin] = cwt(detrend(d_min));
[cmax,fmax] = cwt(detrend(d_max));

% Plot the CWTs
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