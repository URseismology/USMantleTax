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
% Define your clusterCombinationMap
clusterCombinationMap = [1, 1, 4, 4, 4, 4, 1, 4, 4, 4, 1, 3, 3, 3, 2, 3, 2, 3]; %for NVGs
%clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Define the new number of clusters after combination
newNumClusters = length(unique(clusterCombinationMap));

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
        trace_norm(trace_norm > 0) = 0;             %Change here
        modifiedRFs(ii, :) = trace_norm; 

        yvals = trace_norm + ii;
        zeroLine = ii * ones(size(t));

        % Use jbfill to color the positive parts of the trace in blue
        jbfill(t, yvals, zeroLine, 'red', 'none', 0, 1);          %Change here
        plot(t, yvals, 'Color', [0.9, 0.9, 0.9], 'LineWidth', 0.1);
    end

    xlabel('Time (s)');
    ylabel('Station Index');
    xlim([6 30]);
    ylim([1 length(clusterIndices) + 1]); 
    camroll(270); hold off;
   
    %Semblance weighted stacking
    semblance = compute_semblance(modifiedRFs, window_size, step); %Compute semblance for the modified RFs in the current cluster
    % Apply semblance weights to the data
    weightedRFs = modifiedRFs .* repmat(semblance, size(modifiedRFs, 1), 1);
    % Calculate the semblance-weighted mean trace
    meanTrace = nanmean(weightedRFs, 1);
    meanTrace(isnan(meanTrace)) = 0;
    
    %meanTrace = meanTrace - mean(meanTrace); % Centering 
    meanTrace(meanTrace > 0) = 0;      %Change here
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
    fillpart = meanTrace < 0;                  %Change here
    jbfill(t(fillpart), meanTrace(fillpart), zeros(size(meanTrace(fillpart))), 'red', 'none', 0.2, 1);   %Change here

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
%%
clusterNum = 2; % Choose the cluster number

% Retrieve the mean trace for the chosen cluster
meanTraceCluster = meanStacks(clusterNum, :);

% Identify indices of stations in the chosen cluster
clusterIndices = find(newClusters == clusterNum);

% Initialize an array to store modified RFs for visualization
modifiedRFsForVis = zeros(length(clusterIndices), length(RFs(1, :)));

% Apply the same modifications to RFs as in the calculation section
for i = 1:length(clusterIndices)
    stationIndex = clusterIndices(i);
    trace = RFs(stationIndex, :) - mean(RFs(stationIndex, :));
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm; % Apply polarity correction
    trace_norm(trace_norm > 0) = 0;             %change here
    modifiedRFsForVis(i, :) = trace_norm;
end

% Calculate correlations with the mean trace of the cluster
correlations = arrayfun(@(idx) corr(meanTraceCluster', modifiedRFsForVis(idx, :)'), 1:length(clusterIndices));

% Define lower and upper bounds for highlighting correlations
lowerBound = 0.1;
upperBound = max(correlations);

% Normalize correlations for coloring
% Only consider correlations above lowerBound for coloring
normCorrelations = (correlations - lowerBound) / (upperBound - lowerBound);
normCorrelations(normCorrelations < 0) = 0;  % Set values below lowerBound to 0
normCorrelations(normCorrelations > 1) = 1;  % Cap values at 1

% Define the colormap
cmap = customColormap; % Use a standard colormap
colorIndices = round(normCorrelations * (size(cmap, 1) - 1)) + 1;

% Map Projection
figure(300);clf;
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]); %continental US
shading flat;

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Plotting each station in the chosen cluster
for i = 1:length(clusterIndices)
    stationIndex = clusterIndices(i);
    m_line(stationLons(stationIndex), stationLats(stationIndex), 'marker', 'o', 'color', 'k', 'linewi', 1, 'MarkerFaceColor', cmap(colorIndices(i),:), 'MarkerSize', 9);
end

% Coastlines and Grids
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% Adjust colorbar to reflect the correlation range
colormap(customColormap);
caxis([lowerBound, upperBound]);
cb=colorbar('southoutside');
cb.FontSize = 17;

% Title and Labels
title(sprintf('Station Locations for Cluster %d', clusterNum));
xlabel('Longitude');
ylabel('Latitude');

m_text(-127, 23, sprintf('Mean Correlation: %.2f', mean(correlations)), 'FontSize', 15, 'BackgroundColor', 'white');
%print(figure(300),'/Users/evets/Desktop/Earthscope/Figures/N2_stations','-vector','-dpdf','-r0');