addpath(genpath(pwd));
clear;clc;
% Load datasets
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
% Define your clusterCombinationMap
%clusterCombinationMap = [1, 1, 4, 4, 4, 4, 1, 4, 4, 4, 1, 3, 3, 3, 2, 3, 2, 3]; %for NVGs
clusterCombinationMap = [1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4]; %for PVGs

% Create newClusters by mapping original clusters to new clusters
newClusters = arrayfun(@(x) clusterCombinationMap(x), clusters);

% Directly specify colors for the new clusters
newControlPoints = [1, 0, 1;  %magenta                    
                    0, 1, 0  % Green
                    0.6, 0.4, 0.2; % Brown 
                    1, 1, 0];  % Blue
                       
%%
% Create a dendrogram for the entire dataset
figure;
[H, T, outperm] = dendrogram(Z, 0);

% Hide x-axis labels and tick marks
set(gca, 'XTick', [], 'XTickLabel', []);

% Map from original order to new cluster assignments
clusterMapping = newClusters(outperm);

% Color the branches based on the new cluster assignments
idx = 0;
for i = 1:length(H)
    xData = get(H(i), 'XData');
    leafIdx = xData(1);  % Get index of the leaf node

    if i == 1
        idx = 1;
    else
        prevXData = get(H(i-1), 'XData'); % Retrieve XData for the previous line
        if leafIdx ~= prevXData(1)
            idx = idx + 1;
        end
    end

    if idx <= length(clusterMapping)
        newClusterIdx = clusterMapping(idx);
        set(H(i), 'Color', newControlPoints(newClusterIdx, :));
    end
end