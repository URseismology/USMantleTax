addpath(genpath(pwd));
clear; clc;

% Load the datasets
positivesData = readtable('positives_data.txt');
crustThicknessData = readmatrix('US_Crust_GRL_2015_CrustThickness.txt');

% Filter for Cluster 3 stations
cluster3Stations = positivesData(positivesData.Cluster == 3, :);

% Initialize array to store depth differences
depthDifferences = zeros(height(cluster3Stations), 1);

% Loop through each Cluster 3 station
for i = 1:height(cluster3Stations)
    stationLat = cluster3Stations.Latitude(i);
    stationLon = cluster3Stations.Longitude(i);

    % Calculate distances to all crust stations
    distances = sqrt((crustThicknessData(:,1) - stationLat).^2 + ...
                     (crustThicknessData(:,2) - stationLon).^2);

    % Find the index of the closest station
    [~, closestIdx] = min(distances);

    % Calculate depth difference
    depthDifferences(i) = abs(cluster3Stations.Depth(i) - crustThicknessData(closestIdx, 3));
end

% Calculate mean depth difference
meanDepthDifference = mean(depthDifferences);

% Define colors based on the new depth difference condition
colors = zeros(length(depthDifferences), 3); % Initialize color array

% Color criteria
colors(depthDifferences <= 20, :) = repmat([0.5 0.5 0.5], sum(depthDifferences <= 20), 1); % Grey for <= 20
colors(depthDifferences > 20 & depthDifferences <= 40, :) = repmat([0 0 1], sum(depthDifferences > 20 & depthDifferences <= 40), 1); % Blue for > 20 and <= 40
colors(depthDifferences > 40, :) = repmat([1 0 0], sum(depthDifferences > 40), 1); % Red for > 40

% Create scatter plot
figure(1); clf;
scatter(1:length(depthDifferences), depthDifferences, 45, ...
        'filled', 'MarkerEdgeColor', 'k', ...
        'CData', colors); % Set custom colors
xlabel('Station Index');
ylabel('Depth Difference (km)');
title('Depth Difference between Cluster 3 Stations and Closest Crust Thickness Stations');

ax = gca; % Get the current axes handle
ax.FontSize = 20;


% Annotate with mean depth difference
text(1, max(depthDifferences)*0.9, ...
     ['\mu_{depth} = ' num2str(meanDepthDifference, '%.2f') ' km'], ...
     'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'white');

% Draw now to update the figure
drawnow;

% Add annotation with a rectangular border
annotation('rectangle', [0.1, 0.1, 0.8, 0.8], 'LineWidth', 1.5);
print(figure(2),'/Users/evets/Desktop/Earthscope/Figures/S6b','-vector','-dpdf','-r0');
%%
% Section to create histogram of depth differences
figure(2);clf; % Create a new figure for the histogram
histogram(depthDifferences, 'BinWidth', 5); % Adjust 'BinWidth' as needed for better visualization
xlabel('Depth Difference (km)');
ylabel('Frequency');
title('Histogram of Depth Differences for Cluster 3 Stations');
% Get the current axes and increase the font size
ax = gca; 
ax.FontSize = 26; 
