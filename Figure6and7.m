addpath(genpath(pwd));
clear;clc;

% Load seismic receiver functions datasets
unsequencedData = load('RFDataWithCoords_pos_new.mat');
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
% Define the number of clusters
% Calculate pairwise correlation distances and perform clustering
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');

% Determine cluster assignments based on the cutoff distance
cutoffDistance = 0.93893;

clusters = cluster(Z, 'Cutoff', cutoffDistance, 'Criterion', 'distance');
%clusters = cluster(Z, 'Criterion', 'distance','maxclust',7);

% Number of distinct clusters (different colored trees)
NumClusters = numel(unique(clusters));

% Plot dendrogram with a specified cutoff distance for visualization
figure(1);
dendrogram(Z, 0, 'ColorThreshold', cutoffDistance);
title(['Dendrogram with Cutoff Distance = ', num2str(cutoffDistance)]);
xlabel('Sample index');
ylabel('Distance');
set(gca,'XTickLabels','')
hold on;
% Add a horizontal line for the cutoff distance
yline(cutoffDistance, '--r', 'LineWidth', 1.5, 'Label', ['Cutoff = ', num2str(cutoffDistance)]);
hold off;

%%
% %Merge closely related clusters based on time/depth of signal
% clusterCombinationMap = [2, 2, 3, 3, 4, 4, 4, 1]; %For NVG
% newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);
% newNumClusters = length(unique(clusterCombinationMap));

%%
window_size = 10; % window size for semblance computation
step = 1; % step size for semblance computation

% Iterate over each cluster to plot RFs and their semblance-weighted stacks
for clusterNum = 1:NumClusters; %change here
    figure(clusterNum + 1); % Each cluster gets a new figure
    clf; 
    
    % Subplot for individual RFs
    ax1 = subplot(4, 4, [5:4:13, 7:4:15]); hold on;
    clusterIndices = find(clusters == clusterNum); %change here to clusters 
    
    modifiedRFs = zeros(length(clusterIndices), length(RFs(1, :)));

    for ii = 1:length(clusterIndices)
        stationIndex = clusterIndices(ii);
        trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
        trace_norm = trace / max(abs(trace));
        trace_norm = -trace_norm; % Apply polarity correction
        
        % Set all negative amplitudes to zero
        trace_norm(trace_norm > 0) = 0;             %Change here
        modifiedRFs(ii, :) = trace_norm; 

        yvals = trace_norm + ii;
        zeroLine = ii * ones(size(t));

        % Use jbfill to color the positive parts of the trace in blue
        jbfill(t, yvals, zeroLine, 'red', 'none', 0, 1);          %Change here
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);
    end
    
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]); 
    camroll(270); hold off;
    title(['Cluster ', num2str(clusterNum)]);
    xlabel('Time (s)');
    ylabel('Normalized Amplitude + Station Index');
    
    % Semblance-weighted stacking for the current cluster
    semblance = compute_semblance(modifiedRFs, window_size, step); % Compute semblance
    weightedRFs = modifiedRFs .* repmat(semblance, size(modifiedRFs, 1), 1);
    meanTrace = nanmean(weightedRFs, 1);
    meanTrace(isnan(meanTrace)) = 0;

    meanTrace(meanTrace > 0) = 0;      %Change here
    
    % Subplot for mean stack plot
    ax2 = subplot(4, 4, 8:4:16); hold on;
    % Plotting the mean trace 
    plot(t, meanTrace, 'k', 'LineWidth', 1);
    fillpart = meanTrace < 0;                  %Change here
    jbfill(t(fillpart), meanTrace(fillpart), zeros(size(meanTrace(fillpart))), 'red', 'none', 0, 1);   %Change here

    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([6 30]);
    ylim([-1, 1]);
    hold off;
    camroll(270);
end
