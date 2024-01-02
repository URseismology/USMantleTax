% Add the current directory and subdirectories to the path
addpath(genpath(pwd));

%Load datasets
load('negatives_full.mat');  % or negatives_full.mat
unsequencedData = load('RFDataWithCoords_negs_new.mat', 'RFStruct');
sequencedData = load('negatives_full.mat');  % positives dataset

% load('full.mat');  % or negatives_full.mat
% unsequencedData = load('RFDataWithCoords_full_new.mat', 'RFStruct');
% sequencedData = load('full.mat');  %negatives dataset

% Define time vector and data matrix for sequenced and unsequenced data
t = time_vector;
RFs_sequenced = dataset_imageRFs_pos_reordered;

%%
%Plot sequenced data with polarity correction
h1 = figure(1); clf;
hold on;
for ii = 1:size(RFs_sequenced, 1)
    trace = RFs_sequenced(ii, :) - mean(RFs_sequenced(ii, :));
    trace_norm = trace / max(abs(trace));
      
    trace_norm = -trace_norm;

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;
    
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    plot(t, yvals, 'k');
end
xlim([6 30]);
ylim([1 417]);
xlabel('Time (s)');
ylabel('Ordered Station Index');
camroll(270);
%print(h1,'/Users/stevecarr/Desktop/Tolu/wholeUS_filtered_positives_nocrust','-vector','-dpdf','-r0');
%%
%Extract and plot unsequenced data
numTraces = length(unsequencedData.RFStruct);
imageRFs_unsequenced = zeros(numTraces, length(t));
for ii = 1:numTraces
    imageRFs_unsequenced(ii, :) = unsequencedData.RFStruct(ii).Data;
end

h2 = figure(2); clf;
hold on;
for ii = 1:numTraces
    trace_norm = imageRFs_unsequenced(ii, :) / max(abs(imageRFs_unsequenced(ii, :)));
    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    positives = trace_norm > 0;
    negatives = trace_norm < 0;
    
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 0);
    plot(t, yvals, 'k');
end
xlim([6 30]);
ylim([1 numTraces]);
xlabel('Time (s)');
ylabel('Unsequenced Station Index');
camroll(270);
%print(h2,'/Users/stevecarr/Desktop/Tolu/wholeUS_filtered_full_nocrust','-vector','-dpdf','-r0');
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
% Plot locations of sequenced data on a map
lats = [unsequencedData.RFStruct.Lat];
lons = [unsequencedData.RFStruct.Lon];
[~, unsequencedIndices] = ismember({unsequencedData.RFStruct.GridNumber}, {unsequencedData.RFStruct.GridNumber});
[~, sequencedIndices] = ismember({unsequencedData.RFStruct.GridNumber}, cellstr(sequencedData.grid_numbers_reordered));

% Define the latitude and longitude ranges for Contiguous US
uslatlim = [25, 50];
uslonlim = [-125, -67];

% Define the range of indexes you want and plot them on the map
startIndex = 1;
endIndex = 384;
plot_indexes_on_map(startIndex, endIndex, lats, lons, sequencedIndices, uslatlim, uslonlim, customColormap);
%%
% Plot specific station's location and display its indices
gridNumber = 'HRV';
plot_station_location(gridNumber, unsequencedData, sequencedData, uslatlim, uslonlim);
%%
%Make scatter plot of positive and negative depths
% Load positives dataset
load('positives_full.mat');
unsequencedDataPos = load('RFDataWithCoords_pos_new.mat', 'RFStruct');
sequencedDataPos = dataset_imageRFs_pos_reordered;
t = time_vector;

% Create a map for positives max peak times
maxPositivesTimesMap = containers.Map();
for ii = 1:size(sequencedDataPos, 1)
    trace = sequencedDataPos(ii, t >= 8);  % only consider data from 8 seconds upward
    [~, idx] = max(abs(trace));
    startIndex = find(t >= 8, 1);  % find the first index where t >= 6
    maxTime = t(startIndex + idx - 1);
    gridNumber = grid_numbers_reordered(ii, :);
    maxPositivesTimesMap(gridNumber) = maxTime;
end

% Load negatives dataset
load('negatives_full.mat');
unsequencedDataNeg = load('RFDataWithCoords_negs_new.mat', 'RFStruct');
sequencedDataNeg = dataset_imageRFs_pos_reordered;

% Create a map for negatives max peak times
maxNegativesTimesMap = containers.Map();
for ii = 1:size(sequencedDataNeg, 1)
    trace = sequencedDataNeg(ii, :);
    [~, idx] = max(abs(trace));
    maxTime = t(idx);
    gridNumber = grid_numbers_reordered(ii, :);  % Assuming grid_numbers_reordered exists for negatives too
    maxNegativesTimesMap(gridNumber) = maxTime;
end

% Identify common stations
commonGridNumbersPos = keys(maxPositivesTimesMap);
commonGridNumbersNeg = keys(maxNegativesTimesMap);
commonGridNumbers = intersect(commonGridNumbersPos, commonGridNumbersNeg);

% Extract max peak times for common stations
maxPositivesTimesCommon = cell2mat(values(maxPositivesTimesMap, commonGridNumbers));
maxNegativesTimesCommon = cell2mat(values(maxNegativesTimesMap, commonGridNumbers));

% Step 1: Scatter Plot

% Create a scatter plot where points are colored based on the comparison
% between positive and negative peak times.
h4=figure(4);
differences = maxPositivesTimesCommon - maxNegativesTimesCommon;
colors = zeros(length(differences), 3);
colors(differences > 0, :) = repmat([0 0 1], sum(differences > 0), 1);  % Blue for Positives > Negatives
colors(differences <= 0, :) = repmat([1 0 0], sum(differences <= 0), 1);  % Red for Negatives >= Positives
scatter(maxPositivesTimesCommon, maxNegativesTimesCommon, 50, colors, 'filled');
xlabel('PVG Time [s]','FontSize',25);
ylabel('NVG Time [s]','FontSize',25);
title('Detected PVG Time Vs NVG Time','FontSize',25);
xlim([7 33])

% Separate data for reds and blues
reds_x = maxPositivesTimesCommon(differences <= 0);
reds_y = maxNegativesTimesCommon(differences <= 0);

blues_x = maxPositivesTimesCommon(differences > 0);
blues_y = maxNegativesTimesCommon(differences > 0);

% Linear fit for reds
p_reds = polyfit(reds_x, reds_y, 1);
x_fit_reds = linspace(min(reds_x), max(reds_x), 100);
y_fit_reds = polyval(p_reds, x_fit_reds);

% Linear fit for blues
p_blues = polyfit(blues_x, blues_y, 1);
x_fit_blues = linspace(min(blues_x), max(blues_x), 100);
y_fit_blues = polyval(p_blues, x_fit_blues);

hold on;
h_reds = plot(x_fit_reds, y_fit_reds, 'r-', 'LineWidth', 2); % Handle for red line
h_blues = plot(x_fit_blues, y_fit_blues, 'b-', 'LineWidth', 2); % Handle for blue line

legend([h_reds, h_blues], 'NVG>PVG', 'PVG>NVG', 'Location', 'best','FontSize',20); % Use handles to specify which plots to include in the legend
hold off;
%print(h4,'/Users/stevecarr/Desktop/Tolu/scatter_plot','-vector','-dpdf','-r0');

% Step 2: Station Locations with Different Peak Orderings

figure;
usamap('conus');
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none');
hold on;

positiveGreaterThanNegative = differences > 0;

% Extract coordinates for all common stations
allLats = zeros(1, length(commonGridNumbers));
allLons = zeros(1, length(commonGridNumbers));

for ii = 1:length(commonGridNumbers)
    idx = find(strcmp({unsequencedDataPos.RFStruct.GridNumber}, commonGridNumbers{ii}), 1);
    if ~isempty(idx) && isscalar(idx)
        allLats(ii) = unsequencedDataPos.RFStruct(idx).Lat;
        allLons(ii) = unsequencedDataPos.RFStruct(idx).Lon;
    else
        allLats(ii) = NaN;  % Assign NaN if there's no match or multiple matches
        allLons(ii) = NaN;
    end
end

% Use logical indexing to separate the stations
lats_positiveGreater = allLats(positiveGreaterThanNegative);
lons_positiveGreater = allLons(positiveGreaterThanNegative);

lats_negativeGreater = allLats(~positiveGreaterThanNegative);
lons_negativeGreater = allLons(~positiveGreaterThanNegative);

% Plot
geoshow(lats_positiveGreater, lons_positiveGreater, 'DisplayType', 'point', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g','MarkerSize',10);
geoshow(lats_negativeGreater, lons_negativeGreater, 'DisplayType', 'point', 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r','MarkerSize',10);

title('Station Locations: Positives vs Negatives Max Peak Time');

% Step 3: Time Difference Map in station locations

figure();
usamap('conus');
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none');
hold on;

colors_diff = jet(256);  % Example colormap

% Calculate the color limits based on the maximum absolute difference
cLim = max(abs(differences));

% Plot each station with a color indicating the difference
scatterm(allLats, allLons, 80, differences, 'filled', 'MarkerEdgeColor', 'k');

colorbar;
caxis([-cLim, cLim]);  % Set color limits to be symmetric around 0
title('Time Difference between Max Peak Positive and Negative');
colormap(colors_diff);

% Step 4: Time Difference Map with Colormap in finer bins

h6=figure(6);
usamap('conus');

% Remove stations with NaN coordinates
validStations = ~isnan(allLats) & ~isnan(allLons);
validLats = allLats(validStations);
validLons = allLons(validStations);
validDifferences = differences(validStations);

% Define a grid over the US with 0.5-degree bins and extend the grid slightly
[lonGrid, latGrid] = meshgrid((min(validLons)-1):0.8:(max(validLons)+1), ...
                              (min(validLats)-1):0.8:(max(validLats)+1));

% Use griddata for spatial interpolation using 'linear' method
interpolatedDifferences = griddata(validLons, validLats, validDifferences, lonGrid, latGrid, 'linear');

% Load US state boundaries
states = shaperead('usastatelo', 'UseGeoCoords', true);

% Create a binary mask using US state boundaries
mask = false(size(latGrid));

for k = 1:length(states)
    in = inpolygon(lonGrid, latGrid, states(k).Lon, states(k).Lat);
    mask = mask | in;
end

% Mask the interpolated differences
interpolatedDifferences(~mask) = NaN;

% Display the interpolated data as a colormap
pcolorm(latGrid, lonGrid, interpolatedDifferences);
hold on
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none');

colorbar('southoutside');
colormap("jet")
title('Thickness Between PVG and NVG (PVG Time - NVG Time)','FontSize',25);
%print(h6,'/Users/stevecarr/Desktop/Tolu/thickness_map','-vector','-dpdf','-r0');

%Step 5
% Plot Histograms
figure;clf;

% Histogram for Positives
subplot(2,1,1);
histogram(maxPositivesTimesCommon, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Maximum Peak Time (s)');
ylabel('Number of Stations');
title('Distribution of Max Peak Times (Positives)');
grid on;

% Histogram for Negatives
subplot(2,1,2);
histogram(maxNegativesTimesCommon, 'FaceColor', 'r', 'EdgeColor', 'k');
xlabel('Maximum Peak Time (s)');
ylabel('Number of Stations');
title('Distribution of Max Peak Times (Negatives)');
grid on;
%%
%Plot stations based on whether they have a deeper MLD or PVG

% Specify the color of stations you want to plot: 'red' or 'green'
colorToPlot = 'green';  % Change this value to 'green' if you want to plot green stations

% Depending on the specified color, filter the grid numbers
if strcmp(colorToPlot, 'red')
    selectedGridNumbers = commonGridNumbers(differences <= 0);
    plotColor = 'r';
    titleText = 'Locations of Stations where Positive Peak Time <= Negative Peak Time';
elseif strcmp(colorToPlot, 'green')
    selectedGridNumbers = commonGridNumbers(differences > 0);
    plotColor = 'g';
    titleText = 'Locations of Stations where Positive Peak Time > Negative Peak Time';
else
    error('Invalid color specified. Choose "red" or "green".');
end

% Extract corresponding latitudes and longitudes for the selected stations
selectedLats = zeros(size(selectedGridNumbers));
selectedLons = zeros(size(selectedGridNumbers));
for ii = 1:length(selectedGridNumbers)
    idx = find(strcmp({unsequencedDataPos.RFStruct.GridNumber}, selectedGridNumbers{ii}), 1);
    if ~isempty(idx) && isscalar(idx)
        selectedLats(ii) = unsequencedDataPos.RFStruct(idx).Lat;
        selectedLons(ii) = unsequencedDataPos.RFStruct(idx).Lon;
    else
        selectedLats(ii) = NaN;  % Assign NaN if there's no match or multiple matches
        selectedLons(ii) = NaN;
    end
end

% Plot the locations on a map of the US
figure;clf;
usamap('conus');
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'g', 'LineWidth', 1, 'FaceColor', 'none');
hold on;
geoshow(selectedLats, selectedLons, 'DisplayType', 'point', 'Color', plotColor, 'Marker', 'o', 'LineWidth', 1.5);
title(titleText);

%%
% Plot wiggle plots showing the picked maximum amplitude for each trace 
% Set to 1 if positives dataset, 0 if negatives dataset.
isPositiveDataset = 0;  % Change to 0 if negatives dataset

figure(); clf;
hold on;

% Determine the starting index based on the dataset
if isPositiveDataset
    startIndex = find(t >= 8, 1);  % Find the first index where time is 8 seconds or more
else  % Assuming this is the negative dataset
    startIndex = 1;
end

for ii = 1:size(RFs_sequenced, 1)
    trace = RFs_sequenced(ii, :) - mean(RFs_sequenced(ii, :));
    trace_norm = trace / max(abs(trace));
      
    trace_norm = -trace_norm;

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;
    
    % Find the time of the maximum amplitude, starting from the determined index
    [~, maxIndexFromSubset] = max(abs(trace_norm(startIndex:end)));
    maxIndex = maxIndexFromSubset + startIndex - 1;  % Adjusting for the offset introduced by startIndex
    maxTime = t(maxIndex);
    
    % Plot a dashed line at the time of the maximum amplitude
    plot([maxTime, maxTime], [ii, ii + trace_norm(maxIndex)], '--', 'Color', [0, 0, 1], 'LineWidth', 2);
    
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    plot(t, yvals, 'k');
end

xlim([6 30]);
ylim([1 417]);
xlabel('Time (s)');
ylabel('Ordered Station Index');
camroll(270);
%%
%K means Cluster Analysis
numTraces = length(unsequencedDataNeg.RFStruct);
features = zeros(numTraces, 1);  %using max amplitude time

% Identifier for dataset type
datasetType = 'negatives';  % Change to 'positives' or 'negatives' depending on dataset

% K-means Cluster Analysis
numTraces = length(unsequencedDataNeg.RFStruct);

% Create a matrix where each row represents a trace
dataMatrix = zeros(numTraces, length(unsequencedDataNeg.RFStruct(1).Data));

for ii = 1:numTraces
    trace = unsequencedDataNeg.RFStruct(ii).Data;
    
    % If working with positives dataset, only consider data from 8 seconds upward
    if strcmp(datasetType, 'positives')
        trace = trace(t >= 8);
    end

    dataMatrix(ii, :) = trace;
end

numClusters = 3;  % Example number of clusters
[idx, C] = kmeans(dataMatrix, numClusters);

% Custom colors for clustered data
colors_clustered = [
    1 0 0;  % Red
    0 1 0;  % Green
    0 0 1   % Blue
    % 1 1 0  % yellow
];

%Make a map for the clustered data
figure;
usamap('conus');
for ii = 1:numTraces
    lat = unsequencedDataNeg.RFStruct(ii).Lat;
    lon = unsequencedDataNeg.RFStruct(ii).Lon;
    plotm(lat, lon, 'o', 'Color', colors_clustered(idx(ii), :), 'MarkerFaceColor', colors_clustered(idx(ii), :),'MarkerEdgeColor', 'k','MarkerSize', 8);
end
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'g', 'LineWidth', 1, 'FaceColor', 'none');
title(['K-means Clustered data: K =' num2str(numClusters)],'FontSize',20);

% Assuming you've loaded the correct dataset that matches the `datasetType` variable
% Reorder the traces based on the cluster assignments
[~, sortOrder] = sort(idx);
reorderedTraces = sequencedDataNeg(sortOrder, :);

% Plot reordered traces as wiggle plots of K-means cluster
figure();
hold on;

for ii = 1:numTraces
    trace = reorderedTraces(ii, :) - mean(reorderedTraces(ii, :));
    trace_norm = trace / max(abs(trace));
      
    trace_norm = -trace_norm;

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;
    
    % Use jbfill to color the traces based on cluster assignment
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), colors_clustered(idx(sortOrder(ii)), :), 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), colors_clustered(idx(sortOrder(ii)), :), 'none', 0);
    
    plot(t, yvals, 'k');
end

xlim([t(1) t(end)]);
ylim([1 numTraces]);
xlabel('Time (s)');
ylabel('Ordered Station Index (by K-means cluster)');
camroll(270);


%%
% Create a map sequenced data from Sequencer
%Extract latitudes and longitudes from the sequenced data
lats_sequenced = [unsequencedData.RFStruct(sequencedIndices).Lat];
lons_sequenced = [unsequencedData.RFStruct(sequencedIndices).Lon];

figure;
usamap(uslatlim, uslonlim);
states = shaperead('usastatelo', 'UseGeoCoords', true, 'BoundingBox', [uslonlim', uslatlim']);
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'g', 'LineWidth', 1, 'FaceColor', 'none');
hold on;

% Custom colors for sequenced data
color1 = [1, 0, 0]; % Red
color2 = [0, 1, 0]; % Green
color3 = [0, 0, 1]; % Blue
% color4 = [1, 1, 0]; % Yellow

% Plot stations with custom colors
geoshow(lats_sequenced(1:25), lons_sequenced(1:25), 'DisplayType', 'point', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color2);
geoshow(lats_sequenced(26:300), lons_sequenced(26:300), 'DisplayType', 'point', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color1);
geoshow(lats_sequenced(301:end), lons_sequenced(301:end), 'DisplayType', 'point', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color3);

title('Station Locations colored by index range using Sequencer','FontSize',20);
%%




