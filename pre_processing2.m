%% Housekeeping 
close all
clear 
clc

%% Processing for Temperature
% Get all files in the directory
folder = "rems_data";
listing = dir(folder);
% Remove listings that are folders
listing([listing.isdir]) = [];
out_folder = "output";
writematrix(["SOL", "MINIMUM", "MAXIMUM", "MEDIAN", "MEAN", "MODE", "STD"], out_folder + filesep + "temp_stats.csv",'WriteMode', 'overwrite');
writematrix(["SOL", "MINIMUM", "MAXIMUM", "MEDIAN", "MEAN", "MODE", "STD"], out_folder + filesep + "pressure_stats.csv",'WriteMode', 'overwrite');
writematrix(["SOL", "MINIMUM", "MAXIMUM", "MEDIAN", "MEAN", "MODE", "STD"], out_folder + filesep + "humidity_stats.csv",'WriteMode', 'overwrite');
parfor i = 1:length(listing)
    fprintf("Processing file: %s\n", listing(i).name);
    A = readtable(folder + filesep + listing(i).name);
    try
        sol = A.SOL(1);
    catch
        continue;
    end
    temp = A.AMBIENT_TEMP;
    p = A.PRESSURE;
    h = A.LOCAL_RELATIVE_HUMIDITY;
    % Remove NaN values
    temp = temp(~isnan(temp));
    p = p(~isnan(p));
    h = h(~isnan(h));
    temp_data = [sol, min(temp), max(temp), median(temp), mean(temp), mode(temp), std(temp)];
    pres_data = [sol, min(p), max(p), median(p), mean(p), mode(p), std(p)];
    humid_data = [sol, min(h), max(h), median(h), mean(h), mode(h), std(h)];
    if ~any(isnan(temp_data)) && (length(pres_data) == 7)
        writematrix(temp_data, out_folder + filesep + "temp_stats.csv",'WriteMode', 'append');
    end
    if ~any(isnan(pres_data)) && (length(pres_data) == 7)
        writematrix(pres_data, out_folder + filesep + "pressure_stats.csv",'WriteMode', 'append');
    end
    if ~any(isnan(humid_data)) && (length(humid_data) == 7)
        writematrix(humid_data, out_folder + filesep + "humidity_stats.csv",'WriteMode', 'append');
    end
end