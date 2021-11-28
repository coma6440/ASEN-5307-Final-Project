function f = plot_summary(A)
%Remove outliers from data
A = rmoutliers(A);
f = figure('Visible','off');
%Plot the resulting values
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