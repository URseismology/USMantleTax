%Analysis of sequenced Ps seismic receiver functions
addpath(genpath(pwd));
clear;clc;
%%
% Load datasets
load('preferredCentroids.mat', 'preferredCentroids');
unsequencedData = load('RFDataWithCoords_negs_new.mat');
sequencedData = load('negatives_full.mat'); %or change here
stationNames = cellstr(sequencedData.grid_numbers_reordered);

t = sequencedData.time_vector;
RFs_sequenced = sequencedData.dataset_imageRFs_pos_reordered;
RFs_unsequenced = unsequencedData.imageRFs;

%Station coordinates 
stationLats = sequencedData.coordinates_reordered(:, 1); 
stationLons = sequencedData.coordinates_reordered(:, 2);

% Choose the dataset: 'sequenced' or 'unsequenced'
dataChoice = 'sequenced';  

if strcmp(dataChoice, 'unsequenced')
    RFs = RFs_unsequenced;  % Use unsequenced data
else
    RFs = RFs_sequenced;  % Use sequenced data
end
%%
% Creating custom colormap
controlPoints = [0, 0, 0; ...   % black
                 0.5, 0, 0.5; ... % deep purple
                 1, 0.5, 0; ...   % orange
                 1, 1, 0.5];      % light yellow

% Create a linearly spaced vector for interpolation
x = linspace(1, size(controlPoints, 1), 256);  % 256 is the desired number of colors in the colormap

% Interpolate RGB values to create the colormap
customColormap = interp1(1:size(controlPoints, 1), controlPoints, x', 'linear');
%%
% Define the number of clusters
numClusters = 18;                        %change here

% Calculate pairwise correlation distances and perform clustering 
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');
clusters = cluster(Z, 'maxclust', numClusters);

% Sort the traces by cluster assignment
[~, sortOrder] = sort(clusters);

%%
% Plot the clustered receiver functions colored by initial cluster
cmap = lines(numClusters);

% If lines does not provide enough colors, fill the remaining with random colors
if size(cmap, 1) < numClusters
    additionalColors = rand(numClusters - size(cmap, 1), 3);
    cmap = [cmap; additionalColors];
end
controlPoints = cmap;

prevCluster = clusters(sortOrder(1)); % Initialize with the first cluster
clusterStartIndex = 1; % Start index of the current cluster

%%
% Cluster Combination map to combine related clusters
clusterCombinationMap = [1, 1, 1, 4, 4, 4, 1, 4, 2, 4, 1, 3, 3, 3, 2, 4, 2, 4]; %for NVGs
%clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Define the new number of clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Specify colors for the new clusters
newControlPoints = [1, 0, 0;  %Red for cluster 1
                    0.6, 0.4, 0.2; % Brown for cluster 2
                    0, 1, 0;  % Green for cluster 3
                    ];  % Blue for cluster 4

% newControlPoints = [0.6, 0.4, 0.2; % Brown for cluster 1
%                     1, 0, 0;  %Red for cluster 1
%                     0.7, 0.7, 0.7   %Blue for cluster 3
%                     0, 1, 0];  % Green for cluster 3
                    
                       
% Ensure that the number of colors matches the number of new clusters
if size(newControlPoints, 1) < newNumClusters
    additionalColors = rand(newNumClusters - size(newControlPoints, 1), 3);
    newControlPoints = [newControlPoints; additionalColors];
end

% Sort the traces by new cluster assignment
[~, newSortOrder] = sort(newClusters);

%%
% Initialize a cell array to store indices of stations included in the final mean stacks
includedStationsIndices = cell(newNumClusters, 1);

% Iterate over each cluster
for clusterNum = 1:newNumClusters
    figure(clusterNum + 10); clf;

    clusterIndices = find(newClusters == clusterNum);
    numStationsInCluster = length(clusterIndices);
    modifiedRFs = zeros(numStationsInCluster, length(RFs(1, :)));
    includedInMean = false(numStationsInCluster, 1); % To track which traces are included in the final mean
    totalCorrCoeff = 0; % To sum up the correlation coefficients

    % Calculate initial mean trace
    for ii = 1:numStationsInCluster
        stationIndex = clusterIndices(ii);
        trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
        trace_norm = trace / max(abs(trace));
        trace_norm = -trace_norm; % Apply polarity correction
        trace_norm(trace_norm > 0) = 0;      %change here
        modifiedRFs(ii, :) = trace_norm;
    end
    initialMeanTrace = mean(modifiedRFs, 1);

    % Determine which traces correlate well with the initial mean trace
    includedIndices = [];
    for ii = 1:numStationsInCluster
        corrCoeff = corr(initialMeanTrace', modifiedRFs(ii, :)');
        if corrCoeff > 0.5
            includedInMean(ii) = true;
            includedIndices = [includedIndices, clusterIndices(ii)];
            totalCorrCoeff = totalCorrCoeff + corrCoeff;
        end
    end
    includedStationsIndices{clusterNum} = includedIndices;

    % Calculate new mean trace using only the selected traces
    finalMeanTrace = mean(modifiedRFs(includedInMean, :), 1);
    finalMeanStacks(clusterNum, :) = finalMeanTrace;

    % Calculate mean correlation coefficient
    numIncludedTraces = sum(includedInMean);
    meanCorrCoeff = totalCorrCoeff / numIncludedTraces;

    % Subplot for individual RFs
    ax1 = subplot(4, 4, [5:4:13, 7:4:15]); hold on;

    % Plot and fill the traces that contributed to the final mean stack
    for ii = 1:numStationsInCluster
        yvals = modifiedRFs(ii, :) + ii;
        zeroLine = ii * ones(size(t));

        if includedInMean(ii)
            % Fill the area for traces included in the mean stack
            jbfill(t, yvals, zeroLine, 'red', 'none', 0, 1);       %change here
        end

        % Plot the trace on top
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);
    end

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 numStationsInCluster + 1]);
    camroll(270);
    hold off;


    % Subplot for the final mean stack plot
    ax2 = subplot(4, 4, 8:4:16); hold on;
    plot(t, finalMeanTrace, 'k', 'LineWidth', 1);
    fillpart = finalMeanTrace < 0;                  %change here 
    jbfill(t(fillpart), finalMeanTrace(fillpart), zeros(size(finalMeanTrace(fillpart))), 'red', 'none', 0.2, 1); %change here
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([6 30]);
    ylim([-1, 1]);
    hold off;
    camroll(270);

    % Draw rectangles and calculate percentage
    drawnow;
    annotation('rectangle', get(ax1, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
    annotation('rectangle', get(ax2, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
    %percentageInCluster = (numStationsInCluster / totalStations) * 100;

    % Output the mean correlation coefficient and number of stations contributing to the final mean
    fprintf('Cluster %d: Mean Correlation Coefficient = %.2f, %d stations out of %d total stations contributed to the final mean.\n',...
        clusterNum, meanCorrCoeff, numIncludedTraces, numStationsInCluster);
end

%%
% Plot mean stacks on the same figure
ax4 = figure(45); clf; 

% Define the vertical spacing between the plots
verticalSpacing = 0.9;  % Adjust this value as needed for adequate spacing

% Define the order in which clusters should be plotted
plotOrder = [1, 2, 3, 4]; % Change this array to your desired order

% Hold on for plotting multiple stacks in the same subplot
hold on;

% Iterate through the plotOrder array to plot each mean stack
for i = 1:length(plotOrder)
    clusterNum = plotOrder(i);
    meanTrace = finalMeanStacks(clusterNum, :); % Extract the final mean stack for the current cluster
    
    % Offset each mean stack plot vertically
    offsetMeanTrace = meanTrace + (find(plotOrder == clusterNum) - 1) * verticalSpacing;

    % Define the zero line for jbfill with the same offset
    zeroLine = (find(plotOrder == clusterNum) - 1) * verticalSpacing * ones(size(meanTrace));

    % Plot each mean stack with jbfill
    plot(t, offsetMeanTrace, 'k', 'LineWidth', 1.2);
    fillpart = offsetMeanTrace > 0;                     %change here
    jbfill(t, offsetMeanTrace, zeroLine, 'blue', 'none', 0.2, 1);  % change here
    hold on;
end

% Set additional plot properties
xlabel('Time (s)');
ylabel('Amplitude with Offset');
xlim([6 30]); 
ylim([-0.8, (length(plotOrder) - 1) * verticalSpacing + 0.4]); % Adjust Y-axis limits
hold off;
camroll(270);

% Draw the rectangle around the axes
drawnow; % Ensure the figure and its children are fully rendered
axHandle = gca;
annotation('rectangle', get(axHandle, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
%print(figure(45),'/Users/stevecarr/Desktop/Earthscope/Figures/figure8_pvg_newcentroids','-vector','-dpdf','-r0');
%%
%Create a colormap of clustered station on a map of US to investigate patterns
figure(46);
% Define the grid resolution
gridResolution = 0.5; % Adjust as needed
gridLon = min(stationLons):gridResolution:max(stationLons);
gridLat = min(stationLats):gridResolution:max(stationLats);

% Initialize the grid for cluster assignments
clusterGrid = NaN(length(gridLat), length(gridLon)); 

% Define a threshold distance for considering a station (in degrees)
considerationThreshold = 0.9;

% Function to find the indices of stations within the threshold distance
findStationsWithinThreshold = @(x, y, lons, lats, threshold) ...
    find(sqrt((lons - x).^2 + (lats - y).^2) <= threshold);

% Loop over each grid cell
for i = 1:length(gridLat)
    for j = 1:length(gridLon)
        % Find stations within the threshold distance
        nearbyStationsIndices = findStationsWithinThreshold(gridLon(j), gridLat(i), stationLons, stationLats, considerationThreshold);

        % Initialize a variable to store the indices of nearby included stations
        nearbyIncludedIndices = [];

        % Check for nearby included stations in each cluster
        for clusterNum = 1:newNumClusters
            % Extract the indices for the current cluster
            currentClusterIndices = includedStationsIndices{clusterNum};

            % Find the intersection with nearby stations
            includedNearbyIndices = intersect(nearbyStationsIndices, currentClusterIndices);
            
            % Append to the nearby included indices
            nearbyIncludedIndices = [nearbyIncludedIndices; includedNearbyIndices];
        end

        % Logic to decide the cluster of the grid cell based on nearby included stations
        if ~isempty(nearbyIncludedIndices)
            % Find the most common cluster among included nearby stations
            mostCommonCluster = mode(newClusters(nearbyIncludedIndices));
            clusterGrid(i, j) = mostCommonCluster; % Use the cluster number of included stations
        end
    end
end

% Create a custom colormap based on newControlPoints
customColormap = newControlPoints;

% Plot the cluster grid with m_pcolor
m_pcolor(gridLon, gridLat, clusterGrid);
shading flat; % To remove grid lines in m_pcolor plot
colormap(customColormap); % Use the custom colormap initially

% Initialize the m_map projection
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]);
m_coast('line','color','k','linewidth',1);
m_grid('tickdir','out','linewidth',2,'fontsize',10);

% Add Physiographic Provinces
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 0.75);
end

% Add title or other annotations as needed
title('Cluster Assignments with Geographical Features');
%print(figure(46),'/Users/stevecarr/Desktop/Earthscope/Figures/figure8_pvg_newmap','-vector','-dpdf','-r0');
%%
% Map Projection for showing stations by clusters
figure(300);
m_proj('Mercator', 'lon', [-128 -63], 'lat', [25 49.5]); %continental US

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Plotting each station included in the final mean stacks, colored by their cluster
for clusterNum = 1:newNumClusters
    includedClusterIndices = includedStationsIndices{clusterNum}; % Stations included in the final mean for this cluster
    for i = 1:length(includedClusterIndices)
        stationIndex = includedClusterIndices(i);
        m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'color', 'k', 'linewi', 1, 'MarkerFaceColor', customColormap(clusterNum,:), 'MarkerSize', 10);
    end
end

% Coastlines and Grids
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir', 'out', 'linewidth', 2, 'fontsize', 10);
