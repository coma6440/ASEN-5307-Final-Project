%% Housekeeping 
close all
clear 
clc

%% Labels present in data
labels = ["SOL";
          "TIMESTAMP";
          "LMST";
          "LTST";
          "AMBIENT_TEMP";
          "PRESSURE";
          "HORIZONTAL_WIND_SPEED";
          "VERTICAL_WIND_SPEED";
          "WIND_DIRECTION";
          "LOCAL_RELATIVE_HUMIDITY"];
      
%% Options for loading/output folders
folder = "rems_data";
listing = dir(folder);
listing([listing.isdir]) = [];
out_folder = "output";

%% Pre-processing
for i = 1:length(listing)
    A = readmatrix(folder + filesep + listing(i).name, 'OutputType', 'string');
    temp_labels = labels;
    [r,c] = size(A);
    %Remove blank and missing rows
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
