addpath(genpath(pwd));
clear;clc;
% Load datasets
unsequencedData = load('RFDataWithCoords_negs_new.mat');
sequencedData = load('negatives_full.mat'); %or could be negatives_full.mat
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
% Define a spacing variable
spacing = 1;
% Define an offset for extra space before the first and after the last trace
offset = 5;

% Ensure controlPoints has a color for each cluster
% Define colors for each cluster
controlPoints = [1, 0, 0; ...   % red
                 0, 1, 0; ...   % green
                 0, 0, 1; ...   % blue
                 1, 1, 0];      % yellow
%%
% Standardize the receiver function data
RFs_standardized = zscore(RFs);

% Stability Analysis: Evaluating silhouette values for different cluster numbers
E = evalclusters(RFs_standardized, 'kmeans', 'silhouette', 'KList', [1:10]);
stability = E.CriterionValues; % Higher values suggest better stability

% Choose the number of clusters k based on the silhouette values
[~, k] = max(stability);
%%
tic
% Perform k-means clustering with chosen k
[idx, ~] = kmeans(RFs, 4, 'Distance', 'correlation', 'Replicates', 50,'Start','plus');
toc
% Sort the traces by cluster assignment
[~, sortOrder] = sort(idx);

% Plot the clustered receiver functions
figure(1); clf; hold on;
for ii = 1:size(RFs, 1)
    newIndex = sortOrder(ii);
    trace = RFs_sequenced(newIndex, :);  % Get the trace
    trace = RFs_sequenced(newIndex, :) - mean(RFs_sequenced(newIndex, :));
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm;  % Apply polarity correction

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;

    kval = idx(newIndex);
    colGrp = controlPoints(kval, :);  % Get the color for the cluster

    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colGrp, 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), colGrp, 'none', 0);
    plot(t, yvals, 'k');
end

xlim([6 30]);
ylim([1 size(RFs,1) + 1]);

xlabel('Time (s)');
ylabel('Station Index');
title('Clustered Receiver Functions');
camroll(270);
hold off;

%%
% Silhouette Analysis: Calculate silhouette values for each data point
silhouetteValues = silhouette(RFs, idx, 'correlation');

% Calculate the mean silhouette value for each cluster
uniqueClusters = unique(idx);
meanSilhouetteValues = zeros(length(uniqueClusters), 1);

for i = 1:length(uniqueClusters)
    clusterSilhouetteValues = silhouetteValues(idx == uniqueClusters(i));
    meanSilhouetteValues(i) = mean(clusterSilhouetteValues);
    fprintf('Mean silhouette value for cluster %d is %f\n', uniqueClusters(i), meanSilhouetteValues(i));
end

% Plot the silhouette values
figure(2); clf;
silhouette(RFs, idx, 'correlation');
title('Silhouette Values for k-Means Clustering');

% Sensitivity Analysis: Comparing the cluster assignments from initial k-means runto the cluster assignments from multiple runs with different seeds
stabilityResults = zeros(10, 1); % Assess stability for 10 different seeds
for seed = 1:10
    rng(seed);
    [temp_idx, ~] = kmeans(RFs, 4, 'Distance', 'correlation', 'Replicates', 50);
    stabilityResults(seed) = adjustedRandIndex(idx, temp_idx);
end
averageARI = mean(stabilityResults); % Overall stability measure

%%
% Calculate and plot the mean trace for each cluster
figure(3); clf; hold on;
meanTraces = zeros(k, size(RFs, 2));
for i = 1:4
    meanTraces(i, :) = mean(RFs(idx == i, :), 1);
end

verticalOffset = max(abs(meanTraces(:))) * 2; 
for i = 1:4
    offsetMeanTrace = meanTraces(i, :) + (i-1) * verticalOffset;
    plot(t, offsetMeanTrace, 'Color', controlPoints(i, :), 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('Amplitude with Offset');
title('Mean Trace of Each Cluster');
%legend(arrayfun(@(x) ['Cluster ' num2str(x)], 1:k, 'UniformOutput', false));
hold off;
camroll(270)
hold off;

%%
% Plot station locations colored by cluster assignment
% Define clusters to plot and their colors
clustersToPlot = 1:4;
sameColorClusters = [2 3 4]; 
sameColor = [0, 1, 1]; 

% Plot station locations colored by cluster assignment
figure(5); clf; usamap('conus'); hold on;
for i = clustersToPlot % Iterate only over the specified clusters
    clusterIdx = find(idx == i);
    % Check if the current cluster is in the list of clusters to plot with the same color
    if ismember(i, sameColorClusters)
        scatterm(stationLats(clusterIdx), stationLons(clusterIdx), 50, sameColor, 'filled');
    else
        scatterm(stationLats(clusterIdx), stationLons(clusterIdx), 50, controlPoints(i, :), 'filled');
    end
end
hold off;

%%
% Input the station name you want to plot
% stationNameToPlot = 'PKME'; % 
% 
% % Call the function with the necessary parameters
% plotStationLocationByName(stationNameToPlot, stationNames, stationLats, stationLons, idx, controlPoints);
%%
% clusterNumber = 2;  % Replace with the cluster number of interest
% filename = 'stationswithPVG.txt';  % The name of the text file to save the data
% 
% % Call the function with the necessary parameters
% getStationDataByCluster(idx, stationNames, stationLats, stationLons, clusterNumber, filename);

%%
% Load datasets
positivesData = load('positives_full.mat');
negativesData = load('negatives_full.mat');

% Extract relevant data
positiveRFs = positivesData.dataset_imageRFs_pos_reordered;
negativeRFs = negativesData.dataset_imageRFs_pos_reordered;

positiveStationNames = cellstr(positivesData.grid_numbers_reordered);
negativeStationNames = cellstr(negativesData.grid_numbers_reordered);

positiveStationLats = positivesData.coordinates_reordered(:, 1);
positiveStationLons = positivesData.coordinates_reordered(:, 2);

negativeStationLats = negativesData.coordinates_reordered(:, 1);
negativeStationLons = negativesData.coordinates_reordered(:, 2);

% Read the list of stations with both NVG and PVG from the text file
stationsWithBoth = readtable('stationswithPVG.txt', 'ReadVariableNames', false);
stationsWithBothNames = stationsWithBoth.Var1;
%%
% Define colors for stations with both NVG and PVG, and for stations with only NVG
colorNVGandPVG = [0, 1, 0]; % Green for stations with both NVG and PVG
colorOnlyNVG = [1, 0, 0]; % Red for stations with only NVG

% Initialize arrays to hold coordinates
latsNVGandPVG = [];
lonsNVGandPVG = [];
latsOnlyNVG = [];
lonsOnlyNVG = [];

% Aggregate coordinates for stations with both NVG and PVG
for i = 1:length(stationsWithBothNames)
    stationName = stationsWithBothNames{i};
    % Find the station in the negative dataset
    stationIndex = find(strcmp(negativeStationNames, stationName));
    if ~isempty(stationIndex)
        latsNVGandPVG = [latsNVGandPVG; negativeStationLats(stationIndex)];
        lonsNVGandPVG = [lonsNVGandPVG; negativeStationLons(stationIndex)];
    end
end

% Aggregate coordinates for stations with only NVG
for i = 1:length(negativeStationNames)
    stationName = negativeStationNames{i};
    if ~ismember(stationName, stationsWithBothNames)
        latsOnlyNVG = [latsOnlyNVG; negativeStationLats(i)];
        lonsOnlyNVG = [lonsOnlyNVG; negativeStationLons(i)];
    end
end

% Initialize the figure for plotting
figure(6); hold on; grid on;
usamap('conus'); % Set the map to continental US

% Plot stations with both NVG and PVG
scatterm(latsNVGandPVG, lonsNVGandPVG, 50, colorNVGandPVG, 'filled');

% Plot stations with only NVG
scatterm(latsOnlyNVG, lonsOnlyNVG, 50, colorOnlyNVG, 'filled');

%%
% Specify clusters of interest
clustersOfInterest = [2, 3, 4]; % Example clusters

% Initialize arrays for time differences, coordinates, and peak times
timeDifferences = [];
latsForPlotting = [];
lonsForPlotting = [];
peakNVGTimes = [];
peakPVGTimes = [];

% Loop through each station in stationsWithBothNames
for i = 1:length(stationsWithBothNames)
    stationName = stationsWithBothNames{i};
    
    % Find the station in the negative and positive datasets
    negIndex = find(strcmp(negativeStationNames, stationName), 1);
    posIndex = find(strcmp(positiveStationNames, stationName), 1);

    % Check if valid indices were found and if the station is part of the clusters of interest
    if ~isempty(negIndex) && ~isempty(posIndex) && negIndex <= length(idx) && ismember(idx(negIndex), clustersOfInterest)
        % Compute time of maximum amplitude for both datasets
        [~, negMaxIdx] = max(abs(negativeRFs(negIndex, :)));
        [~, posMaxIdx] = max(abs(positiveRFs(posIndex, :)));

        % Compute the time difference (NVG time - PVG time)
        negMaxTime = negativesData.time_vector(negMaxIdx);
        posMaxTime = positivesData.time_vector(posMaxIdx);
        timeDifference = posMaxTime - negMaxTime;

        % Store time difference, peak times, and coordinates for plotting
        timeDifferences = [timeDifferences; timeDifference];
        peakNVGTimes = [peakNVGTimes; negMaxTime];
        peakPVGTimes = [peakPVGTimes; posMaxTime];
        latsForPlotting = [latsForPlotting; negativeStationLats(negIndex)];
        lonsForPlotting = [lonsForPlotting; negativeStationLons(negIndex)];
    end
end

% Plot the time differences on a map
figure(7); clf;hold on; 
usamap('conus'); % Set the map to continental US
scatterm(latsForPlotting, lonsForPlotting, 50, timeDifferences, 'filled');
colormap(jet); % Set the colormap to 'jet' for better visualization
colorbar; % Add a colorbar to indicate the time differences
title('Time Difference between NVG and PVG (NVG time - PVG time)');
hold off;

% Plot histograms of peak NVG and PVG times
figure(8);clf;
subplot(2,1,1);
histogram(peakNVGTimes, 'BinWidth', 0.5); % Adjust BinWidth as needed
title('Histogram of Peak NVG Times');
xlabel('Time (s)');
ylabel('Frequency');

subplot(2,1,2);
histogram(peakPVGTimes, 'BinWidth', 0.5); % Adjust BinWidth as needed
title('Histogram of Peak PVG Times');
xlabel('Time (s)');
ylabel('Frequency');
%%
% % Initialize figure and subplots
% figure(10);clf;
% subplot(2,1,1); % NVG subplot
% hold on;
% title('Wiggle Plots of NVGs with Picked Peak Time');
% xlabel('Time (s)');
% ylabel('Station Index');
% 
% % Plot NVG wiggle plots with peak times using jbfill
% for i = 1:length(stationsWithBothNames)
%     stationName = stationsWithBothNames{i};
%     negIndex = find(strcmp(negativeStationNames, stationName), 1);
%     if ~isempty(negIndex) && ismember(idx(negIndex), clustersOfInterest)
%         trace = negativeRFs(negIndex, :);
%         trace_norm = trace / max(abs(trace));
%         trace_norm = -trace_norm;  % Apply polarity correction
%         yvals = trace_norm + i; % Adjust the vertical position for each trace
%         negatives = trace_norm < 0;
%         positives = trace_norm > 0;
% 
%         % Use jbfill to plot
%         jbfill(negativesData.time_vector(negatives), yvals(negatives), i * ones(size(negatives)), 'r', 'none', 0);
%         jbfill(negativesData.time_vector(positives), yvals(positives), i * ones(size(positives)), 'r', 'none', 0);
% 
%         % Mark peak time
%         [~, negMaxIdx] = max(abs(trace));
%         plot(negativesData.time_vector(negMaxIdx), yvals(negMaxIdx), 'k.', 'MarkerSize', 8, 'LineWidth', 2);
%     end
%     xlim([5 30]);
%     camroll(90)
% end
% hold off;
% 
% subplot(2,1,2); % PVG subplot
% hold on;
% title('Wiggle Plots of PVGs with Picked Peak Time');
% xlabel('Time (s)');
% ylabel('Station Index');
% 
% subplot(2,1,2); % PVG subplot
% hold on;
% title('Wiggle Plots of PVGs with Picked Peak Time');
% xlabel('Time (s)');
% ylabel('Station Index');
% 
% % Plot PVG wiggle plots with peak times using jbfill
% for i = 1:length(stationsWithBothNames)
%     stationName = stationsWithBothNames{i};
%     posIndex = find(strcmp(positiveStationNames, stationName), 1);
%     if ~isempty(posIndex) && ismember(idx(posIndex), clustersOfInterest)
%         trace = positiveRFs(posIndex, :);
%         trace_norm = trace / max(abs(trace));
%         trace_norm = -trace_norm;  % Apply polarity correction
%         yvals = trace_norm + i; % Adjust the vertical position for each trace
%         negatives = trace_norm < 0;
%         positives = trace_norm > 0;
% 
%         % Use jbfill to plot
%         jbfill(positivesData.time_vector(negatives), yvals(negatives), i * ones(size(negatives)), 'b', 'none', 0);
%         jbfill(positivesData.time_vector(positives), yvals(positives), i * ones(size(positives)), 'b', 'none', 0);
% 
%         % Mark peak time
%         [~, posMaxIdx] = max(abs(trace));
%         plot(positivesData.time_vector(posMaxIdx), yvals(posMaxIdx), 'k.', 'MarkerSize', 8, 'LineWidth', 2);
%     end
%     xlim([5 30]);
%     camroll(90)
% end
% hold off;


