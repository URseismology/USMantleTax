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
%clusterCombinationMap = [1, 1, 4, 4, 4, 4, 1, 4, 4, 4, 1, 3, 3, 3, 2, 3, 2, 3]; %for NVGs
clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Define the new number of clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

% Directly specify colors for the new clusters
newControlPoints = [1, 0, 1;  %magenta                    
                    0, 1, 0  % Green
                    0.6, 0.4, 0.2; % Brown 
                    1, 1, 0];  % Blue
                       
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
    %meanTrace = pws2(modifiedRFs, p); %phase weighted stacking
    %meanTrace = f_pws2(modifiedRFs, p, window_size, step); %phase and semblance weighted stacking

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
    %fprintf('Cluster %d has %.2f%% out of %d total stations.\n', clusterNum, percentageInCluster, totalStations);
end
%print(figure(14),'/Users/evets/Desktop/Earthscope/Figures/figure6d','-vector','-dpdf','-r0');
