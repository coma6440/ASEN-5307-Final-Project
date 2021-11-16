%% Housekeeping 
close all
clear 
clc

%% Loading in data
labels = ["SOL";
          "TIMESTAMP";
          "LMST";
          "LTST";
          "AMBIENT_TEMP";
          "PRESSURE";
          "HORIZONTAL_WIND_SPEED";
          "VERTICAL_WIND_SPEED";
          "WIND_DIRECTION";
          "WS_CONFIDENCE_LEVEL";
          "PRESSURE_UNCERTAINTY";
          "VOLUME_MIXING_RATIO";
          "LOCAL_RELATIVE_HUMIDITY"];
% Get all files in the directory
folder = "rems_data";
listing = dir(folder);
% Remove listings that are folders
listing([listing.isdir]) = [];
out_folder = "output";
%Got to 1813
for i = 51:length(listing)
    A = readmatrix(folder + filesep + listing(i).name, 'OutputType', 'string');
    temp_labels = labels;
    [r,c] = size(A);
    A = deblank(A);
    A(A == "") = missing;
    for j = 5:c
        idx = ~ismissing(A(:,j));
        if any(idx)
            writematrix([str2double(A(idx,1)),...
                         str2double(A(idx,2)),...
                         A(idx,3),...
                         A(idx,4),...
                         str2double(A(idx,j))],...
                         out_folder + filesep + labels(j) + ".csv",...
                         'WriteMode', 'append');
        end
    end
end
