% Add the current directory and subdirectories to the path
addpath(genpath(pwd));

% Load datasets
myData = readtable('NVG_depth_data2.csv');
Abt_data = readtable('Abt_NVG_2010.txt');
Liu_data = readtable('Liu_NVG_Combined.csv');
HUA_data = readtable('HUA_MLD.txt');

% Define the threshold for proximity (1 degree)
proximity_threshold = 1; % degrees

% Initialize arrays for storing scatter plot data and markers
scatterData = [];
scatterMarkers = {};

% Define markers for each study
markers = {'^', 'p', 's'}; % Triangle for Abt, Pentagon for HUA, Square for Liu/Shear (shallowMLD or LAB), Diamond for Liu/Shear (deepMLD)

% List of studies and their datasets
studies = {'Abt', 'Liu', 'HUA'};
study_datasets = {Abt_data, Liu_data, HUA_data};

% Process each study
for i = 1:length(studies)
    study_data = study_datasets{i};

    % Set column names based on the dataset
    if strcmp(studies{i}, 'Liu')
        lat_col = 'Latitude';
        long_col = 'Longitude';
        depth_col = 'Depth';
        interp_col = 'interp';
    elseif strcmp(studies{i}, 'HUA')
        lat_col = 'latitude';
        long_col = 'longitude';
        depth_col = 'NVG_depth_km_';
    else % For Abt
        lat_col = 'Lat';
        long_col = 'Long';
        depth_col = 'Depth_km_';
    end

    for j = 1:height(study_data)
        studyLat = study_data{j, lat_col};
        studyLong = study_data{j, long_col};
        studyDepth = study_data{j, depth_col};

        % Special handling for Liu data
        if strcmp(studies{i}, 'Liu')
            studyInterp = study_data{j, interp_col};
            if strcmp(studyInterp, 'deepMLD')
                myDataFiltered = myData(myData.Cluster == 2 | myData.Cluster == 3, :);
                marker = 'h'; % Diamond for deepMLD
            else % shallowMLD or LAB
                myDataFiltered = myData(myData.Cluster == 1, :);
                marker = 'p'; % Square for shallowMLD or LAB
            end
        else
            myDataFiltered = myData;
            marker = markers{i};
        end

        distances = sqrt((myDataFiltered.Latitude - studyLat).^2 + (myDataFiltered.Longitude - studyLong).^2);
        [minDistance, idx] = min(distances);
        if minDistance <= proximity_threshold
            myDepth = myDataFiltered.Depth(idx);

            %if myDepth >= 60 && myDepth <= 150
                depthDifference = abs(myDepth - studyDepth);
                scatterData = [scatterData; myDepth, studyDepth, depthDifference];
                scatterMarkers{end+1} = marker;
            %end
        end
    end
end

% Calculate the full range of depth data for x-axis
fullRangeMin = min(scatterData(:, 1));
fullRangeMax = max(scatterData(:, 1));

% Set the y-axis limit
yAxisLimit = 160;
%%
% Initialize variables for colored data
coloredScatterData = [];

% Plotting all comparisons on one figure
fig = figure(1); clf;
hold on;

% Plot scatter data with conditional coloring
for k = 1:size(scatterData, 1)
    if scatterData(k, 3) > 25
        % Plot data with depth difference greater than 30 km in light grey
        scatter(scatterData(k, 1), scatterData(k, 2), 70, [0.8 0.8 0.8], 'filled', 'Marker', scatterMarkers{k}, 'MarkerEdgeColor', 'k');
    else
        % Plot data with depth difference less than or equal to 30 km with depth difference coloring
        scatter(scatterData(k, 1), scatterData(k, 2), 70, scatterData(k, 3), 'filled', 'Marker', scatterMarkers{k}, 'MarkerEdgeColor', 'k');
        coloredScatterData = [coloredScatterData; scatterData(k, :)];
    end
end

% Add a one-to-one line
plot([fullRangeMin, yAxisLimit], [fullRangeMin, yAxisLimit], 'r--', 'LineWidth', 1.5);

% Customize plot
colormap('jet');
caxis([min(coloredScatterData(:, 3)), max(coloredScatterData(:, 3))]);
colorbar('southoutside');
xlabel('My Depth (km)');
ylabel('Study Depth (km)');
title('Depth Comparisons with Previous Studies');
set(gca, 'FontSize', 20);
xlim([fullRangeMin, fullRangeMax]);
ylim([fullRangeMin, yAxisLimit]);

% Add the legend
%legend([hAbt, hLiuShallow, hLiuDeep, hHua], {'Abt', 'Liu (shallowMLD/LAB)', 'Liu (deepMLD)', 'HUA'}, 'Location', 'northeast');

% Calculate R-squared and mean difference for colored plots
rSquared = fitlm(coloredScatterData(:, 1), coloredScatterData(:, 2)).Rsquared.Ordinary;
meanDiff = mean(coloredScatterData(:, 3));
textStr = sprintf('R^2: %.2f\nMean Absolute Deviation: %.2f km', rSquared, meanDiff);
text(mean([fullRangeMin, fullRangeMax]), yAxisLimit - 10, textStr, 'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

hold off;

% Rectangle annotation
% axHandle = gca;
% drawnow;
% annotation('rectangle', get(axHandle, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
