function [min_fit, max_fit, f] = plot_fit(t, data, min_fit, max_fit, yaxis_str)
%Evaluate the fit over t
min_fit = feval(min_fit,t);
max_fit = feval(max_fit,t);

%Plot the data and resulting fit
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