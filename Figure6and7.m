addpath(genpath(pwd));
clear;clc;
% Load datasets
load('preferredCentroids.mat', 'preferredCentroids');
unsequencedData = load('RFDataWithCoords_negs_new.mat');
sequencedData = load('positives_full.mat'); %or change here
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
numClusters = 11;        %change here

% Generate a colormap with distinct colors
cmap = lines(numClusters);

% If lines does not provide enough colors, fill the remaining with random colors
if size(cmap, 1) < numClusters
    additionalColors = rand(numClusters - size(cmap, 1), 3);
    cmap = [cmap; additionalColors];
end
controlPoints = cmap;
%%
% Calculate pairwise correlation distances and perform clustering
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');
clusters = cluster(Z, 'maxclust', numClusters);

% Sort the traces by cluster assignment
[~, sortOrder] = sort(clusters);

%%
% Plot the clustered receiver functions
figure(1); clf; hold on;
prevCluster = clusters(sortOrder(1)); % Initialize with the first cluster
clusterStartIndex = 1; % Start index of the current cluster

for ii = 1:size(RFs, 1)
    newIndex = sortOrder(ii);
    currentCluster = clusters(newIndex);

    % Check for a change in cluster
    if currentCluster ~= prevCluster
        % Output the index range of the previous cluster
        fprintf('Cluster %d: Station Index Range %d to %d\n', prevCluster, clusterStartIndex, ii-1);
        
        % Update for the next cluster
        prevCluster = currentCluster;
        clusterStartIndex = ii;
    end

    trace = RFs(newIndex, :);  % Get the trace
    trace = RFs(newIndex, :) - mean(RFs(newIndex, :));
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm;  % Apply polarity correction

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;

    kval = clusters(newIndex);
    if kval == -1
        kval = numClusters + 1;  % Assuming the last color is for outliers
    elseif kval > numClusters
        kval = mod(kval, numClusters) + 1;  % Ensure kval doesn't exceed the number of defined colors
    end

    colGrp = controlPoints(kval, :);

    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colGrp, 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), colGrp, 'none', 0);
    plot(t, yvals, 'k');
end

% Output for the last cluster
fprintf('Cluster %d: Station Index Range %d to %d\n', prevCluster, clusterStartIndex, size(RFs, 1));

xlim([6 30]);
ylim([1 size(RFs,1) + 1]);
xlabel('Time (s)');
ylabel('Station Index');
camroll(270);
hold off;

%%
% Define your clusterCombinationMap
%clusterCombinationMap = [1, 1, 1, 4, 4, 4, 1, 4, 2, 4, 1, 3, 3, 3, 2, 4, 2, 4]; %for NVGs
clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Define the new number of clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Directly specify colors for the new clusters
newControlPoints = [1, 0, 0;   % Red
                    1, 1, 0;   % Yellow
                    0, 1, 0;   % Green
                    0, 0, 1];  % Blue

% Ensure that the number of colors matches the number of new clusters
if size(newControlPoints, 1) < newNumClusters
    additionalColors = rand(newNumClusters - size(newControlPoints, 1), 3);
    newControlPoints = [newControlPoints; additionalColors];
end

% Sort the traces by new cluster assignment
[~, newSortOrder] = sort(newClusters);


% Plot the newly clustered receiver functions
figure(2); clf; hold on;
for ii = 1:size(RFs, 1)
    newIndex = newSortOrder(ii);
    trace = RFs_sequenced(newIndex, :);  % Get the trace
    trace = trace - mean(trace);
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm;  % Apply polarity correction

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;

    kval = newClusters(newIndex);
    colGrp = newControlPoints(kval, :);

    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colGrp, 'none',0,1);
    jbfill(t(positives), yvals(positives), zeroLine(positives), colGrp, 'none',0,1);
    plot(t, yvals, 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.5);
end

xlim([6 30]);
ylim([1 size(RFs,1) + 1]);
xlabel('Time (s)');
ylabel('Station Index');
camroll(270);
hold off;

%%
% Define the number of new clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Define a matrix to store mean traces for each cluster
meanStacks = zeros(newNumClusters, length(t));  % Assuming 't' is the time vector

for clusterNum = 1:newNumClusters
    figure(clusterNum + 10); clf; % Create a new figure for each new cluster

    % Subplot for individual RFs
    ax1 = subplot(4, 4, [5:4:13, 7:4:15]); hold on;
    clusterIndices = find(newClusters == clusterNum); % Find indices in the new cluster

    % Initialize an array to store modified RFs for each cluster
    modifiedRFs = zeros(length(clusterIndices), length(RFs(1, :)));

    for ii = 1:length(clusterIndices)
        stationIndex = clusterIndices(ii);
        trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
        trace_norm = trace / max(abs(trace));
        trace_norm = -trace_norm; % Apply polarity correction
        
        % Set all negative amplitudes to zero
        trace_norm(trace_norm < 0) = 0;             %Change here
        modifiedRFs(ii, :) = trace_norm; 

        yvals = trace_norm + ii;
        zeroLine = ii * ones(size(t));

        % Use jbfill to color the positive parts of the trace in blue
        jbfill(t, yvals, zeroLine, 'blue', 'none', 0, 1);     %Change here
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);
    end

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]); 
    camroll(270); hold off;
   
    % Calculate mean trace from modified RFs
    meanTrace = mean(modifiedRFs, 1);
    meanTrace = meanTrace - mean(meanTrace); % Centering 
    meanTrace(meanTrace < 0) = 0;      %Change here
    meanStacks(clusterNum, :) = meanTrace;  % Store the mean trace

    % Subplot for mean stack plot
    ax2 = subplot(4, 4, 8:4:16); hold on;
    % Plotting the mean trace 
    plot(t, meanTrace, 'k', 'LineWidth', 1);
    fillpart = meanTrace > 0;                  %Change here
    jbfill(t(fillpart), meanTrace(fillpart), zeros(size(meanTrace(fillpart))), 'blue', 'none', 0.2, 1);   %Change here

    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([6 30]);
    ylim([-1, 1]);
    hold off;
    camroll(270);
    
    % Draw rectangles after the subplots have been created
    % Use annotation to draw rectangles
    drawnow; % Ensure the figure and its children are fully rendered
    annotation('rectangle', get(ax1, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
    annotation('rectangle', get(ax2, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
end
%print(figure(14),'/Users/evets/Desktop/Earthscope/Figures/figure7d','-vector','-dpdf','-r0');

%%
%figure 8
% Define the clusters you want to highlight
clustersToHighlight = [];
highlightColor = [1, 1, 1];  % White color for highlighted clusters

% Map Projection for Figure 44
figure(44); clf;
 % Subplot for station locations
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]);
m_coast('line','color','k','linewidth',1);
m_grid('tickdir','out','linewidth',2,'fontsize',10);

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Plot station locations colored by new cluster assignment
for i = 1:newNumClusters  % Iterate over the new clusters
    clusterIdx = find(newClusters == i);  % Find indices of stations in the new cluster

    for j = 1:length(clusterIdx)
        stationIndex = clusterIdx(j);

        % Check if the current cluster is one of those to be highlighted
        if ismember(i, clustersToHighlight)
            m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', highlightColor, 'MarkerSize', 12, 'LineStyle', 'none');
        else
            m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', newControlPoints(i, :), 'MarkerSize', 12, 'LineStyle', 'none');
        end
    end
end
%%
% Define the figure for plotting mean stacks
ax4=figure(45); clf; 

% Define the vertical spacing between the plots
verticalSpacing = 0.5;  % Adjust this value as needed for adequate spacing

% Define the order in which clusters should be plotted
plotOrder = [2, 1, 4, 3]; % Change this array to your desired order

% Hold on for plotting multiple stacks in the same subplot
hold on;

% Iterate through the plotOrder array to plot each mean stack
for i = 1:length(plotOrder)
    clusterNum = plotOrder(i);
    meanTrace = meanStacks(clusterNum, :); % Extract the mean stack for the current cluster
    
    % Offset each mean stack plot vertically
    % The order of mean stack is determined by their position in plotOrder
    offsetMeanTrace = meanTrace + (find(plotOrder == clusterNum) - 1) * verticalSpacing;

    % Define the zero line for jbfill with the same offset
    zeroLine = (find(plotOrder == clusterNum) - 1) * verticalSpacing * ones(size(meanTrace));

    % Plot each mean stack with jbfill
    plot(t, offsetMeanTrace, 'k', 'LineWidth', 0.8);
    fillpart = offsetMeanTrace > 0;
    jbfill(t, offsetMeanTrace, zeroLine, 'blue', 'none', 0.2, 1); 
    hold on
end

% Set additional plot properties
xlabel('Time (s)');
ylabel('Amplitude with Offset');
xlim([6 30]); 
ylim([-0.4, (length(plotOrder) - 1) * verticalSpacing + 0.4]); % Adjust Y-axis limits based on the number of clusters and spacing
hold off;
camroll(270)

% Retrieve the axes handle after plotting
axHandle = gca;

% Draw the rectangle around the axes
drawnow; % Ensure the figure and its children are fully rendered
annotation('rectangle', get(axHandle, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
%%
% Define the cluster number to visualize
clusterNum = 1; % Change this to your desired cluster number

timeRangeStart = 12; % change here
timeRangeEnd = 18;   

% Find indices of the time range
timeRangeIndices = find(t >= timeRangeStart & t <= timeRangeEnd);

% Initialize an array to store the time of maximum amplitude for each station
% The length should be the total number of stations, not just the clusterIndices
maxAmplitudeTimes = zeros(size(RFs, 1), 1); % Adjusted to the total number of stations

% Subplot for individual RFs
figure(100); clf; hold on;
clusterIndices = find(newClusters == clusterNum); % Find indices in the cluster

assert(max(clusterIndices) <= size(RFs, 1), 'Cluster indices exceed the number of stations.');

for ii = 1:length(clusterIndices)
    stationIndex = clusterIndices(ii);
    trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm; % Apply polarity correction

    % Set all negative amplitudes to zero
    trace_norm(trace_norm < 0) = 0;     %change here

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));

    % Use jbfill to color the positive parts of the trace in blue
    jbfill(t, yvals, zeroLine, 'blue', 'none', 0, 1);      %change here
    plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);

    % Find and mark the maximum positive amplitude within the time range
    positiveAmplitudes = trace_norm(timeRangeIndices);
    positiveAmplitudes(positiveAmplitudes <= 0) = NaN; % Ignore non-positive values %change here <=
    [~, maxIdx] = max(positiveAmplitudes); %change here to max/min
    maxTime = t(timeRangeIndices(maxIdx));
    maxAmplitudeTimes(stationIndex) = maxTime; % Store the time of maximum amplitude
    plot([maxTime, maxTime], [ii-0.5, ii+0.5], 'k-', 'LineWidth', 1); % Black dashed line
end

xlabel('Time (s)');
ylabel('Station Index');
xlim([6 30]);
ylim([1 length(clusterIndices) + 1]);
camroll(270); hold off;

% Map Projection
figure(200);
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]); %continental US
shading flat;

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Plotting each station in the cluster
cmap = customColormap; % Colormap
for i = 1:length(clusterIndices)
    stationIndex = clusterIndices(i);
    % Check if the time is within the defined range
    if maxAmplitudeTimes(stationIndex) >= timeRangeStart && maxAmplitudeTimes(stationIndex) <= timeRangeEnd
        % Normalize the time to colormap index
        colorIdx = round(interp1(linspace(timeRangeStart, timeRangeEnd, size(cmap,1)), 1:size(cmap,1), maxAmplitudeTimes(stationIndex), 'linear', 'extrap'));
        colorIdx = max(1, min(size(cmap,1), colorIdx)); % Ensure index is within bounds
    else
        % Assign a default color (e.g., gray) if the time is outside the range
        colorIdx = size(cmap,1) + 1; % Index for the default color
        cmap = [cmap; 0.5, 0.5, 0.5]; % Append the default color to the colormap
    end
    m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'color', 'k', 'linewi', 1, 'MarkerFaceColor', cmap(colorIdx,:), 'MarkerSize', 12);
end

% Coastlines and Grids
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% Add colorbar
colormap(cmap);
caxis([timeRangeStart, timeRangeEnd]);
colorbar('southoutside');

% Title and Labels
title(sprintf('Station Locations for Cluster %d', clusterNum));
xlabel('Longitude');
ylabel('Latitude');
%%
% Visualize station locations for a chosen cluster
clusterNum = 2; % Change this to your desired cluster number

% Calculate mean trace and correlations (assuming you have the RFs and newClusters arrays)
clusterIndices = find(newClusters == clusterNum);
modifiedRFs = zeros(length(clusterIndices), length(RFs(1, :)));
for i = 1:length(clusterIndices)
    stationIndex = clusterIndices(i);
    trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
    trace(trace < 0) = 0; % Set negative amplitudes to zero
    modifiedRFs(i, :) = trace;
end
meanTrace = mean(modifiedRFs, 1);

% Calculate variances and correlations
variances = arrayfun(@(i) mean((modifiedRFs(i, :) - meanTrace).^2), 1:length(clusterIndices));
correlations = arrayfun(@(i) corr(meanTrace', modifiedRFs(i, :)'), 1:length(clusterIndices));
meanVariance = mean(variances);
meanCorrelation = mean(correlations);

% Map Projection
figure(300);
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]); %continental US
shading flat;

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Define the colormap
cmap = customColormap; % Use a standard colormap
% Normalize correlations for coloring
normCorrelations = (correlations - min(correlations)) / (max(correlations) - min(correlations));

colorIndices = round(normCorrelations * (size(cmap, 1) - 1)) + 1;

% Plotting each station in the cluster
for i = 1:length(clusterIndices)
    stationIndex = clusterIndices(i);
    m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'color', 'k', 'linewi', 1, 'MarkerFaceColor', cmap(colorIndices(i),:), 'MarkerSize', 12);
end

% Coastlines and Grids
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% Add colorbar
colormap(customColormap);
caxis([min(correlations), max(correlations)]);
colorbar('southoutside');

% Title and Labels
title(sprintf('Station Locations for Cluster %d', clusterNum));
xlabel('Longitude');
ylabel('Latitude');

% Add mean variance and mean correlation as text
m_text(-127, 23, sprintf('$\\bar{r}: %.2f$', meanCorrelation), 'FontSize', 15, 'BackgroundColor', 'white', 'Interpreter', 'latex');
%print(figure(300),'/Users/evets/Desktop/Earthscope/Figures/negatives_3','-vector','-dpdf','-r0');

%%

