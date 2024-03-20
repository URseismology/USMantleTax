addpath(genpath(pwd));
addpath('./m_map')
clear;clc;

% Load datasets
sequencedData = load('positives_full.mat'); %or change here
stationNames = cellstr(sequencedData.grid_numbers_reordered);

t = sequencedData.time_vector;
RFs_sequenced = sequencedData.dataset_imageRFs_pos_reordered;

%Station coordinates 
stationLats = sequencedData.coordinates_reordered(:, 1); 
stationLons = sequencedData.coordinates_reordered(:, 2);

RFs = RFs_sequenced;  % Use sequenced data

%%
% Define the number of clusters
% Calculate pairwise correlation distances and perform clustering
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');

% Determine cluster assignments based on the cutoff distance
 cutoffDistance = 0.94;     

clusters = cluster(Z, 'Cutoff', cutoffDistance, 'Criterion', 'distance');
%%
%Merge closely related clusters based on time/depth of signal
%clusterCombinationMap = [2, 2, 3, 3, 4, 4, 4, 1]; %For NVG
clusterCombinationMap = [4, 2, 4, 3, 3, 1, 2, 1]; %For PVG

newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);
newNumClusters = length(unique(clusterCombinationMap));

% Directly specify colors for the new clusters
% newControlPoints = [1, 0, 0;     % Red
%                     1, 1, 0;     % Yellow
%                     0.35, 0.16, 0.14; % Syrup (Dark Brown)
%                     0.6, 0, 0.8]; %Violet

newControlPoints = [0, 0, 0.5;         % Navy Blue
                    0.53, 0.81, 0.98;  % Sky Blue                   
                    0.94, 0.94, 0.94;     % light grey
                    0, 0.5, 0.5];      % Teal Blue
                    
% Ensure that the number of colors matches the number of new clusters
if size(newControlPoints, 1) < newNumClusters
    additionalColors = rand(newNumClusters - size(newControlPoints, 1), 3);
    newControlPoints = [newControlPoints; additionalColors];
end

%%
%Parameters for semblance weighted stacking
window_size = 10; %semblance window size
step = 1;  %step size for semblance computation

% Define the number of new clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Calculate the total number of stations
totalStations = size(RFs, 1); 

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
        jbfill(t, yvals, zeroLine, 'blue', 'none', 0, 1);          %Change here
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);
    end

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]); 
    camroll(270); hold off;
   
    %Use phase and semblance weighted stacking instead of simple mean
    %meanTrace = mean(modifiedRFs, 1); %simple mean stack

    %Semblance weighted stacking
    semblance = compute_semblance(modifiedRFs, window_size, step); %Compute semblance for the modified RFs in the current cluster
    % Apply semblance weights to the data
    weightedRFs = modifiedRFs .* repmat(semblance, size(modifiedRFs, 1), 1);
    % Calculate the semblance-weighted mean trace
    meanTrace = nanmean(weightedRFs, 1);
    meanTrace(isnan(meanTrace)) = 0;
    
    %meanTrace = meanTrace - mean(meanTrace); % Centering 
    meanTrace(meanTrace < 0) = 0;                       %Change here
    meanStacks(clusterNum, :) = meanTrace;  % Store the mean trace

    % Initialize an array to store correlation coefficients
    correlationCoefficients = zeros(1, length(clusterIndices));

    % Loop through each trace in the cluster
    for ii = 1:length(clusterIndices)
        % Calculate the correlation coefficient between each trace and the mean trace
        correlationCoefficients(ii) = corr(meanTrace', modifiedRFs(ii, :)');
    end

    % Calculate the mean correlation coefficient for the cluster
    meanCorrelation = nanmean(correlationCoefficients);

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
    drawnow; % Ensure the figure and its children are fully rendered
    annotation('rectangle', get(ax1, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
    annotation('rectangle', get(ax2, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);

    % Calculate the percentage of stations in this cluster
    numStationsInCluster = length(clusterIndices);
    percentageInCluster = (numStationsInCluster / totalStations) * 100;

     % Calculate correlation of each trace with mean stack
    corr_coeffs = NaN(size(modifiedRFs,1),1); 
    for i = 1:size(modifiedRFs,1)
        corr_matrix = corrcoef(modifiedRFs(i,:),meanTrace);
        corr_coeffs(i) = corr_matrix(1,2); 
    end

     % Output the mean correlation coefficient for the current cluster
    fprintf('Mean correlation coefficient for cluster %d: %.4f\n', clusterNum, meanCorrelation);

    % Output the result to the command window
    fprintf('Cluster %d has %.2f%% out of %d total stations.\n', clusterNum, percentageInCluster, totalStations);
end

%%
% Load the data
data = readtable('positives_data.txt');

stationLats = data.Latitude; 
stationLons = data.Longitude;
clusters = data.Cluster;

% Define the grid resolution
gridResolution = 0.5; % Adjust as needed
gridLon = min(stationLons):gridResolution:max(stationLons);
gridLat = min(stationLats):gridResolution:50;

% Initialize the grid for cluster assignments
clusterGrid = NaN(length(gridLat), length(gridLon)); 

% Define a threshold distance for considering a station (in degrees)
considerationThreshold = 0.8;

% Function to find the indices of stations within the threshold distance
findStationsWithinThreshold = @(x, y, lons, lats, threshold) ...
    find(sqrt((lons - x).^2 + (lats - y).^2) <= threshold);

% Loop over each grid cell
for i = 1:length(gridLat)
    for j = 1:length(gridLon)
        % Find stations within the threshold distance
        nearbyStationsIndices = findStationsWithinThreshold(gridLon(j), gridLat(i), stationLons, stationLats, considerationThreshold);

        % Logic to decide the cluster of the grid cell based on nearby stations
        if ~isempty(nearbyStationsIndices)
            % Find the most common cluster among nearby stations
            mostCommonCluster = mode(clusters(nearbyStationsIndices));
            clusterGrid(i, j) = mostCommonCluster; % Directly use the cluster number
        end
    end
end

% Initialize base map
figure(46); clf;
m_proj('Mercator','lon',[-126 -63],'lat',[25 49.5]); %continental US
m_coast('patch', [0.78 0.78 0.78], 'edgecolor', 'none'); % Set coast to grey to blend with background
m_grid('tickdir','out','linewidth',2,'fontsize',10);

hold on;
% Plot the cluster grid with m_pcolor
m_pcolor(gridLon, gridLat, clusterGrid);

% Read US geological provinces
S = shaperead('physio.shp');

% List of provinces to plot
provincesToPlot = {'COASTAL PLAIN','PACIFIC BORDER'};

% Loop through the shapes and plot only the specified provinces
for k = 1:length(S)
    if any(strcmp(S(k).PROVINCE, provincesToPlot))
        m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1.5);
    end
end

% % states
M=m_shaperead('ne_50m_admin_1_states_provinces');
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color', [0.6 0.6 0.6],'linewi',0.05);
end

% Read provinces boundary coordinates from the CSV file
Basin_and_range = readtable('Basin_and_range.csv');
Appalatians = readtable('Appalatians.csv');
Snake_river_plain = readtable('Snake_river_plain.csv');
Colorado_plateau = readtable('Colorado_plateau.csv');
Rocky_mountains = readtable('Rocky_mountains.csv');

m_line(Basin_and_range.Longitude, Basin_and_range.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Appalatians.Longitude, Appalatians.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Snake_river_plain.Longitude, Snake_river_plain.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Colorado_plateau.Longitude, Colorado_plateau.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Rocky_mountains.Longitude, Rocky_mountains.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% Create a custom colormap based on newControlPoints
customColormap = newControlPoints;
%shading flat; % To remove grid lines in m_pcolor plot
colormap(customColormap); % Use the custom colormap initially
%%
% Define the figure for plotting mean stacks
ax4=figure(45); clf; 

% Define the vertical spacing between the plots
verticalSpacing = 0.5;  % Adjust this value as needed for adequate spacing

% Define the order in which clusters should be plotted
%plotOrder = [2, 1, 4, 3];
plotOrder = [1, 2, 3, 4];% Change here

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
    plot(t, offsetMeanTrace, 'k', 'LineWidth', 1.2);
    fillpart = offsetMeanTrace > 0;             %change here
    jbfill(t, offsetMeanTrace, zeroLine, 'blue', 'none', 0.2, 1);  %change here
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