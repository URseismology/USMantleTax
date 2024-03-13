addpath(genpath(pwd));
clear;clc;

% Load seismic receiver functions datasets
unsequencedData = load('RFDataWithCoords_pos_new.mat');
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
% Define the number of clusters
numClusters = 11;        %test change here
%%
% Calculate pairwise correlation distances and perform clustering
D = pdist(RFs, 'correlation');
Z = linkage(D, 'average');
clusters = cluster(Z, 'maxclust', numClusters);

% Sort the traces by cluster assignment
[~, sortOrder] = sort(clusters);

%%
% Define your clusterCombinationMap
clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; 

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Define the new number of clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Sort the traces by new cluster assignment
[~, newSortOrder] = sort(newClusters);

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

end
