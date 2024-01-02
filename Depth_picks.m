addpath(genpath(pwd));
clear;clc;

% Load datasets
unsequencedData = load('RFDataWithCoords_negs_new.mat');
sequencedData = load('negatives_full.mat');         % change here
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
% Define the number of clusters
numClusters = 18;        %change here

%%
% Calculate pairwise correlation distances and perform clustering
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');
clusters = cluster(Z, 'maxclust', numClusters);

% Sort the traces by cluster assignment
[~, sortOrder] = sort(clusters);

%%
% Define your clusterCombinationMap
clusterCombinationMap = [1, 1, 4, 4, 4, 4, 1, 4, 4, 4, 1, 3, 3, 3, 2, 3, 2, 3]; %for NVGs
%clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);
                
%%
%Section plots receiver functions belonging to each cluster
% Define the time ranges for each cluster
%timeRanges = {[11, 18], [8, 12], [6, 8], [24, 30]}; % For positives data
timeRanges = {[6, 10], [10, 13], [12, 20], [19, 30]}; % For negatives data


% Define the number of new clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Initialize an array to store time of maximum amplitude
timeOfMaxAmplitudeArray = zeros(size(RFs, 1), 1);

for clusterNum = 1:newNumClusters
    figure(clusterNum + 10); clf; % Create a new figure for each new cluster
    hold on
    clusterIndices = find(newClusters == clusterNum); % Find indices in the new cluster

    % Time range for the current cluster
    currentRange = timeRanges{clusterNum};
    timeIndices = find(t >= currentRange(1) & t <= currentRange(2));

    % Initialize an array to store modified RFs for each cluster
    modifiedRFs = zeros(length(clusterIndices), length(RFs(1, :)));

    % Initialize an array to store non-zero times of max amplitude
    nonZeroTimes = [];

    for ii = 1:length(clusterIndices)
        stationIndex = clusterIndices(ii);
        trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
        trace_norm = trace / max(abs(trace));
        trace_norm = -trace_norm; % Apply polarity correction

        trace_norm(trace_norm > 0) = 0; %change here
        modifiedRFs(ii, :) = trace_norm;

        yvals = trace_norm + ii;
        zeroLine = ii * ones(size(t));

        % Use jbfill to color the positive parts of the trace in blue
        jbfill(t, yvals, zeroLine, 'red', 'none', 0, 1);     %change here
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);

        % Find and mark the highest positive amplitude within the time range
        [maxAmplitude, maxIndex] = min(trace_norm(timeIndices));      %chnage here
        if maxAmplitude < 0          %change here
            timeOfMaxAmplitude = t(timeIndices(maxIndex));
            timeOfMaxAmplitudeArray(clusterIndices(ii)) = timeOfMaxAmplitude;
            nonZeroTimes = [nonZeroTimes, timeOfMaxAmplitude]; % Collect non-zero times
            plot(t(timeIndices(maxIndex)), ii + maxAmplitude, 'k.', 'MarkerSize', 15, 'LineWidth', 2);
        end
    end

     % Calculate the mean of non-zero times of max amplitude for the cluster
    meanTimeOfMaxAmplitude = mean(nonZeroTimes);

    % Assign the mean time to zero values in the cluster
    zeroAmplitudeIndices = clusterIndices(timeOfMaxAmplitudeArray(clusterIndices) == 0);
    timeOfMaxAmplitudeArray(zeroAmplitudeIndices) = meanTimeOfMaxAmplitude;

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]);
    camroll(270); hold off;
    
    ax = gca;
    % Draw rectangles after the subplots have been created
    drawnow; % Ensure the figure and its children are fully rendered
    annotation('rectangle', get(ax, 'Position'), 'EdgeColor', 'k', 'LineWidth', 1);
end
print(figure(14),'/Users/evets/Desktop/Earthscope/Figures/S4d','-vector','-dpdf','-r0');
%%
% Load Schultz-Pelkum dataset
schultzData = readtable('Schultz_Pelkum_vpvs.csv');

% Calculate depths for each station and store in a matrix
depths = zeros(length(RFs), 1); % Initialize array for depths

for i = 1:size(RFs,1)
    stationLat = stationLats(i);
    stationLon = stationLons(i);
    
    % Calculate distances to all points in the Schultz-Pelkum dataset
    distances = sqrt((schultzData.latitude - stationLat).^2 + (schultzData.longitude - stationLon).^2);
    [~, idx] = min(distances);
    
    vs = schultzData.vsv(idx);
    vp = schultzData.vp(idx);
    
    timeOfMaxAmplitude = timeOfMaxAmplitudeArray(i);
    
    % Calculate depth if timeOfMaxAmplitude is not zero (i.e., it was set)
    if timeOfMaxAmplitude > 0
        depth = timeOfMaxAmplitude / (1/vs - 1/vp);
        depths(i) = depth;
    end
end

% Corrections: set depths greater than x to y
cluster1Indices = find(newClusters == 1);
depths(cluster1Indices(depths(cluster1Indices) > 100)) = 99;

cluster2Indices = find(newClusters == 2);
depths(cluster2Indices(depths(cluster2Indices) < 99)) = 99;

cluster3Indices = find(newClusters == 3);
depths(cluster3Indices(depths(cluster3Indices) < 130)) = 130;

cluster3Indices = find(newClusters == 3);
depths(cluster3Indices(depths(cluster3Indices) > 200)) = 199;

cluster4Indices = find(newClusters == 4);
depths(cluster4Indices(depths(cluster4Indices) > 300)) = 300;

% Plot depth histograms for each cluster
for clusterNum = 1:newNumClusters
    figure(clusterNum + 20); % New figure for depth histogram
    clusterIndices = find(newClusters == clusterNum);
    histogram(depths(clusterIndices),'BinWidth',5);
    xlabel('Depth (km)');
    ylabel('Frequency');
    title(['Depth Histogram for Cluster ', num2str(clusterNum)]);
end
%print(figure(21),'/Users/evets/Desktop/Earthscope/Figures/N1_hist','-vector','-dpdf','-r0');
