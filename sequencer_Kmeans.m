addpath(genpath(pwd));

% Load datasets
unsequencedData = load('RFDataWithCoords_negs_new.mat');
%unsequencedData = load('RFDataWithCoords_full_new.mat');
sequencedData = load('negatives_full.mat'); 

% Define time vector and data matrix for sequenced or unsequenced data
t = sequencedData.time_vector;
RFs_sequenced = sequencedData.dataset_imageRFs_pos_reordered;
RFs_unsequenced = unsequencedData.imageRFs;

% Extract station coordinates
stationLats = sequencedData.coordinates_reordered(:, 1); % For sequenced data
stationLons = sequencedData.coordinates_reordered(:, 2);

% Define a spacing variable
spacing = 1;
% Define an offset for extra space before the first and after the last trace
offset = 5;

%%
% Ensure controlPoints has a color for each cluster
% Define colors for each cluster
controlPoints = [1, 0, 0; ...   % red
                 0, 1, 0; ...   % green
                 0, 0, 1; ...   % blue
                 1, 1, 0;
                 0.5 0.5 0.5];      % yellow
%%
% Standardize the receiver function data
RFs_standardized = zscore(RFs_unsequenced);

tic
% Determine the optimal number of clusters using Gap statistic
%optimalK = evalclusters(RFs_standardized, 'kmeans', 'Gap', 'KList', [1:10]);

% Check Davies-Bouldin index for different numbers of clusters
dbIndex = zeros(1, 10);
for k = 1:10
    [idx, ~] = kmeans(RFs_standardized, k, 'Distance', 'correlation');
    dbIndex(k) = daviesbouldin(RFs_standardized, idx);
end
[minDb, bestK_Db] = min(dbIndex);

% Stability Analysis (Consensus Clustering)
E = evalclusters(RFs_standardized, 'kmeans', 'silhouette', 'KList', [1:20]);
stability = E.CriterionValues; % Higher values suggest better stability

toc

% Choose the number of clusters
%numClusters = optimalK.OptimalK;

%%
% Set a random seed for reproducibility
%rng(123);

% Perform k-means clustering
k = 2;
[idx, ~] = kmeans(RFs_unsequenced, k, 'Distance', 'correlation', 'Replicates', 10);

% Sort the traces by cluster assignment
[sortedIdx, sortOrder] = sort(idx);

% Plot the wiggle plots with clusters
h3 = figure(3); clf;
hold on;
for ii = 1:size(RFs_unsequenced, 1)
    newIndex = sortOrder(ii);
    %trace = RFs_sequenced(ii, :) - mean(RFs_sequenced(ii, :));
    trace = RFs_unsequenced(newIndex, :);
    trace_norm = trace / max(abs(trace));
    %trace_norm = -trace_norm;  % Apply polarity correction

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;

    kval = idx(newIndex);
    colGrp = controlPoints(kval, :);  % Get the color for the cluster

    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colGrp, 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    plot(t, yvals, 'k');
end
xlim([6 30]);
ylim([1 size(RFs_unsequenced,1) + 1]);

xlabel('Time (s)');
ylabel('Station Index');
title('Clustered Receiver Functions');
camroll(270);
hold off;

%%
% Calculate the silhouette coefficient for the clustering
silhouetteValues = silhouette(RFs_unsequenced, idx, 'correlation');

% Plot the silhouette values
figure;
silhouette(RFs_unsequenced, idx, 'correlation');
title('Silhouette Values for k-Means Clustering');

% Perform sensitivity analysis
stabilityResults = zeros(10, 1); % Assuming we want to check stability for 10 different seeds
for seed = 1:10
    rng(seed);
    [temp_idx, ~] = kmeans(RFs_unsequenced, k, 'Distance', 'correlation', 'Replicates', 10);
    % Compare temp_idx with idx to check the stability using the Adjusted Rand Index
    stabilityResults(seed) = adjustedRandIndex(idx, temp_idx);
end

% Now you can take the average of the stabilityResults to get an overall stability measure
averageARI = mean(stabilityResults);

% Calculate and plot the mean trace for each cluster
meanTraces = zeros(k, size(RFs_unsequenced, 2));
for i = 1:k
    meanTraces(i, :) = mean(RFs_unsequenced(idx == i, :), 1);
end

% Plot the mean traces with vertical offsets
figure;
hold on;
verticalOffset = max(abs(meanTraces(:))) * 2; % Adjust this value as needed
for i = 1:k
    offsetMeanTrace = meanTraces(i, :) + (i-1) * verticalOffset;
    plot(t, offsetMeanTrace, 'Color', controlPoints(i, :), 'LineWidth', 2);
end
xlabel('Time (s)');
ylabel('Amplitude with Offset');
title('Mean Trace of Each Cluster');
legend(arrayfun(@(x) ['Cluster ' num2str(x)], 1:k, 'UniformOutput', false));
hold off;
camroll(270)
%%
% Plot the station locations colored by their cluster assignment
figure;
usamap('conus');
hold on;
for i = 1:k
    clusterIdx = find(idx == i);
    scatterm(stationLats(clusterIdx), stationLons(clusterIdx), 50, controlPoints(i, :), 'filled');
end
legend(arrayfun(@(x) ['Cluster ' num2str(x)], 1:k, 'UniformOutput', false));
title('Station Locations Colored by Cluster');
hold off;
