addpath(genpath(pwd));
clear;clc;
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
% Define the number of clusters
numClusters = 18;        %change here

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
clusterCombinationMap = [1, 1, 1, 4, 4, 4, 1, 4, 2, 4, 1, 3, 3, 3, 2, 4, 2, 4]; %for NVGs
%clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

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

    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colGrp, 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), colGrp, 'none', 0);
    plot(t, yvals, 'k');
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
        trace_norm(trace_norm > 0) = 0;             %Change here
        modifiedRFs(ii, :) = trace_norm; 

        yvals = trace_norm + ii;
        zeroLine = ii * ones(size(t));

        % Use jbfill to color the positive parts of the trace in blue
        jbfill(t, yvals, zeroLine, 'red', 'k', 0, 0.5);     %Change here
        plot(t, yvals, 'k', 'LineWidth', 1); 
    end

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]); 
    camroll(270); hold off;
   
    % Calculate mean trace from modified RFs
    meanTrace = mean(modifiedRFs, 1);
    meanTrace = meanTrace - mean(meanTrace); % Centering 
    meanTrace(meanTrace > 0) = 0;      %Change here

    % Subplot for mean stack plot
    ax2 = subplot(4, 4, 8:4:16); hold on;
    % Plotting the mean trace 
    plot(t, meanTrace, 'k', 'LineWidth', 1);
    fillpart = meanTrace < 0;                  %Change here
    jbfill(t(fillpart), meanTrace(fillpart), zeros(size(meanTrace(fillpart))), 'red', 'none', 1, 0.5);   %Change here

    % Calculate variances and correlations using the modified RFs
    variances = arrayfun(@(i) mean((modifiedRFs(i, :) - meanTrace).^2), 1:length(clusterIndices));
    correlations = arrayfun(@(i) corr(meanTrace', modifiedRFs(i, :)'), 1:length(clusterIndices));
    meanVariance = mean(variances);
    meanCorrelation = mean(correlations);

    % Add variance and correlation information as text
    text(max(t)*0.7, max(meanTrace) + 0.1, sprintf('$\\bar{\\sigma}^2: %.2f$', meanVariance), 'FontSize', 15, 'BackgroundColor', 'white', 'Interpreter', 'latex');
text(max(t)*0.5, max(meanTrace), sprintf('$\\bar{r}: %.2f$', meanCorrelation), 'FontSize', 15, 'BackgroundColor', 'white', 'Interpreter', 'latex');


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

    % Subplot for histogram of correlations
    % subplot(4, 4, 1:3); 
    % histogram(correlations, 'BinWidth', 0.05); % Adjust BinWidth as needed
    % xlim([-1, 1]); % Correlations range from -1 to 1
    % xlabel('Correlation');
    % ylabel('Frequency');
    % title('Correlation Histogram');
end
%print(figure(14),'/Users/evets/Desktop/Earthscope/Figures/figure6d','-vector','-dpdf','-r0');
%%
% Define the clusters you want to highlight
clustersToHighlight = [];
highlightColor = [1, 1, 1];  % Black color for highlighted clusters

% Plot station locations colored by new cluster assignment
figure(4); clf; usamap('conus'); hold on;

for i = 1:newNumClusters % Iterate over the new clusters
    clusterIdx = find(newClusters == i); % Find indices of stations in the new cluster

    % Check if the current cluster is one of those to be highlighted
    if ismember(i, clustersToHighlight)
        scatterm(stationLats(clusterIdx), stationLons(clusterIdx), 50, highlightColor, 'filled');
    else
        scatterm(stationLats(clusterIdx), stationLons(clusterIdx), 50, newControlPoints(i, :), 'filled');
    end
end

hold off;
