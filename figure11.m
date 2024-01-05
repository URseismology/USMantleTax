% Add the current directory and subdirectories to the path
addpath(genpath(pwd));

% Load datasets
myNVGdata = readtable('NVG_depth_data2.csv');
Abt_data = readtable('Abt_NVG_2010.txt');
Liu_data = readtable('Liu_NVG_Combined.csv');
Hua_data = readtable('HUA_MLD.txt');
Kreuger_data = readtable('Kreuger_UMD.txt');
Hua_PVG_data = readtable('Hua_PVG.txt');
myPVGdata = readtable('positives_data.txt');

% Initialize arrays for storing scatter plot data and markers
scatterData = [];
scatterMarkers = [];
scatterColors = [];

% Define the threshold for proximity (1 degree)
proximity_threshold = 1; % degrees

% Process each study
for study = {'Abt', 'Liu', 'Hua', 'Kreuger', 'Hua_PVG'}
    studyName = study{1};
    study_data = eval([studyName, '_data']);

     % Set column names based on the dataset
    if strcmp(studyName, 'Kreuger')
        lat_col = 'Lat';
        long_col = 'Long';
        depth_col = 'Depth';
        interp_col = 'Interp';
    elseif strcmp(studyName, 'Hua_PVG')
        lat_col = 'latitude';
        long_col = 'longitude';
        depth_col = 'PVG_150_depth_km';
        % No interp column for Hua_PVG
    elseif strcmp(studyName, 'Liu')
        lat_col = 'Latitude';
        long_col = 'Longitude';
        depth_col = 'Depth';
        interp_col = 'interp';
    elseif strcmp(studyName, 'Hua')
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

        % Set studyInterp for Liu and Kreuger, leave it empty for others
        if strcmp(studyName, 'Kreuger') || strcmp(studyName, 'Liu')
            studyInterp = study_data{j, interp_col};
        else
            studyInterp = '';
        end

        % Determine which cluster to compare with based on depth and study
         if strcmp(studyName, 'Liu')
            if strcmp(studyInterp, 'deepMLD')
                myDataFiltered = myNVGdata(myNVGdata.Cluster == 2, :);
                currentCluster = 2;
            else
                myDataFiltered = myNVGdata(myNVGdata.Cluster == 1, :);
                currentCluster = 1;
            end
        elseif strcmp(studyName, 'Hua_PVG') || (strcmp(studyName, 'Kreuger') && strcmp(studyInterp, 'PVG'))
            if studyDepth > 200
                myDataFiltered = myPVGdata(myPVGdata.Cluster == 4, :);
                currentCluster = 4;
            else
                myDataFiltered = myPVGdata(myPVGdata.Cluster == 1, :);
                currentCluster = 1;
            end
        else
            if studyDepth <= 100
                myDataFiltered = myNVGdata(myNVGdata.Cluster == 1, :);
                currentCluster = 1;
            elseif studyDepth > 100 && studyDepth <= 150
                myDataFiltered = myNVGdata(myNVGdata.Cluster == 2, :);
                currentCluster = 2;
            else
                myDataFiltered = myNVGdata(myNVGdata.Cluster == 3 | myNVGdata.Cluster == 4, :);
                currentCluster = 3;
            end
         end
         
        % Calculate distance and compare depths
        distances = sqrt((myDataFiltered.Latitude - studyLat).^2 + (myDataFiltered.Longitude - studyLong).^2);
        [minDistance, idx] = min(distances);

        if minDistance > proximity_threshold && currentCluster < 4           
            % Check the next cluster if no close station is found
            nextCluster = determineNextCluster(currentCluster);
            myDataFiltered = myNVGdata(myNVGdata.Cluster == nextCluster, :);
            distances = sqrt((myDataFiltered.Latitude - studyLat).^2 + (myDataFiltered.Longitude - studyLong).^2);
            [minDistance, idx] = min(distances);
        end

        if minDistance <= proximity_threshold
            myDepth = myDataFiltered.Depth(idx);
            depthDifference = abs(myDepth - studyDepth);
            scatterData = [scatterData; myDepth, studyDepth, depthDifference];

            % Determine the marker and color based on the study
            scatterMarker = determineMarker(studyName, studyInterp);
            scatterMarkers{end+1} = scatterMarker;
            if depthDifference > 30
                scatterColor = [0.7,0.7,0.7]; % Grey color for next cluster
            else
                scatterColor = determineColor(studyName, studyInterp); % Original color
            end
            scatterColors{end+1} = scatterColor;
        end
    end
end

% Calculate the full range of depth data for x-axis
fullRangeMin = min(scatterData(:, 1));
fullRangeMax = max(scatterData(:, 1));
%%
% Plotting all comparisons on one figure
fig = figure(1); clf;
hold on;

% Initialize variables for R-squared calculation
validScatterData = [];

% Plot scatter data and collect valid data for R-squared calculation
for k = 1:length(scatterData)
    if scatterData(k, 3) <= 30  % Check if the depth difference is <= 30
        validScatterData = [validScatterData; scatterData(k, 1:2)]; % Collect only depth data
    end
    if ~isequal(scatterColors{k}, [0.5, 0.5, 0.5])  % Check if the color is not gray
        scatter(scatterData(k, 1), scatterData(k, 2), 70, scatterColors{k}, 'filled', 'Marker', scatterMarkers{k}, 'MarkerEdgeColor', 'k');
    else
        scatter(scatterData(k, 1), scatterData(k, 2), 70, [0.5, 0.5, 0.5], 'filled', 'Marker', scatterMarkers{k}, 'MarkerEdgeColor', 'k');
    end
end

% Add a one-to-one line starting from 50
axisMin = 50;
fullRangeMax = max(max(scatterData(:, 1:2)));
plot([axisMin, fullRangeMax], [axisMin, fullRangeMax], 'k--', 'LineWidth', 1.5);

% Customize plot
xlabel('My Depth (km)');
ylabel('Study Depth (km)');
title('Depth Comparisons with Previous Studies');
set(gca, 'FontSize', 20);
xlim([axisMin, fullRangeMax]);
ylim([axisMin, fullRangeMax]);

% Calculate R-squared and mean difference for plots with depth difference <= 30
rSquared = fitlm(validScatterData(:, 1), validScatterData(:, 2)).Rsquared.Ordinary;
meanDiff = mean(scatterData(scatterData(:, 3) <= 30, 3)); % Mean difference for depth difference <= 30
textStr = sprintf('R^2: %.2f\nMean Absolute Deviation: %.2f km', rSquared, meanDiff);
text(fullRangeMax * 0.3, fullRangeMax * 0.9, textStr, 'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');

hold off;

% Rectangle annotation
axHandle = gca;
drawnow;
annotation('rectangle', get(axHandle, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
%print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/figure9x','-vector','-dpdf','-r0');
%%
% Histogram for depth differences
fig2 = figure(2); clf;
histogram(scatterData(:, 3)); % Assuming the 3rd column is the depth difference

% Customize histogram
xlabel('Depth Difference (km)');
ylabel('Frequency');
title('Histogram of Depth Differences');
set(gca, 'FontSize', 20);

% Optional: set specific bins for the histogram
% binEdges = 0:10:max(scatterData(:, 3)); % Change the bin size as needed
% histogram(scatterData(:, 3), 'BinEdges', binEdges);
%print(figure(2),'/Users/evets/Desktop/Earthscope/Figures/figure9xx','-vector','-dpdf','-r0');