addpath(genpath(pwd));
clear; clc;

% Load the data
NVG_data = readtable('NVG_depth_data2.csv');
positives_data = readtable('positives_data.txt');

% Initialize arrays for storing class 1 and class 2 station data
class1Stations = [];
class1Depths = [];
class2Stations = [];
class2Shapes = []; % 'o' for circle, 's' for square
class2Thicknesses = [];

% Classify stations
for i = 1:height(NVG_data)
    station = NVG_data(i, :);
    station_name = station.StationName{1};
    station_cluster = station.Cluster;
    station_depth = station.Depth;

    idx_in_positives = find(strcmp(positives_data.StationName, station_name));

    % Class 1 logic
    if isempty(idx_in_positives) && (station_cluster == 1 || station_cluster == 2)
        class1Stations = [class1Stations; station];
        class1Depths = [class1Depths; station_depth];
    elseif ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
        if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 3)
            class1Stations = [class1Stations; station];
            class1Depths = [class1Depths; station_depth];
        end
    end


    % Class 2 logic
    if ~isempty(idx_in_positives) && (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 1 | cluster_in_positives == 2)
        for j = 1:length(idx_in_positives)
            depth_in_positives = positives_data.Depth(idx_in_positives(j));
            thickness = abs(station_depth - depth_in_positives);
            class2Stations = [class2Stations; station];
            class2Thicknesses = [class2Thicknesses; thickness];
            if station_depth < depth_in_positives
                class2Shapes = [class2Shapes; 'o']; % Circle for lower depth
            else
                class2Shapes = [class2Shapes; 's']; % Square for higher depth
            end
        end
    end
end
%%
% % Creating custom colormap
controlPoints = [0, 0, 0; ...   % black
                 0.5, 0, 0.5; ... % deep purple
                 1, 0.5, 0; ...   % orange
                 1, 1, 0.5];      % light yellow

% Create a linearly spaced vector for interpolation
x = linspace(1, size(controlPoints, 1), 256);  % 256 is the desired number of colors in the colormap

% Interpolate RGB values to create the colormap
customColormap1 = interp1(1:size(controlPoints, 1), controlPoints, x', 'linear');

% Custom colormap from mud brown to deep blue
controlPoints = [
    0.40, 0.26, 0.13; ... % dirt brown
    0.65, 0.50, 0.39; ... % yellowish brown
    0.68, 0.85, 0.90; ... % light green
    0.53, 0.81, 0.98; ... % sky blue
    0.00, 0.00, 0.55  ... % deep blue
];

% Create a linearly spaced vector for interpolation
x = linspace(1, size(controlPoints, 1), 256);  % 256 is the desired number of colors in the colormap

% Interpolate RGB values to create the colormap
customColormap2 = interp1(1:size(controlPoints, 1), controlPoints, x, 'linear');

%stationCmap = customColormap;
%stationCmap = flipud(gray(256));
stationCmap = flipud(customColormap1);
%%
% Load seismic velocity data
velocityData = readtable('Schultz_Pelkum_vpvs.csv');

% Plot Class 1 Stations Colored by Depth
figure(1); clf;

% Create main axes for plotting
mainAx = axes('Position', [0.1, 0.3, 0.6, 0.6]);

% Plot seismic velocities using m_scatter
m_scatter(velocityData.longitude, velocityData.latitude, 250, velocityData.vsv, 'Marker', '.');
caxis([min(velocityData.vsv), max(velocityData.vsv)]); % Adjust as needed
colormap(mainAx,customColormap2);

cb1 = colorbar(mainAx,'eastoutside');
cb1.Position = [0.75, 0.45, 0.03, 0.5];
cb1.FontSize = 13;
freezeColors(mainAx);  % Freeze the colors of the plot
cbfreeze(cb1); % Freeze the colors of the first colorbar
ylabel(cb1, 'Velocity');

m_proj('Mercator', 'lon', [-128 -63], 'lat', [25 49.5]);
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir', 'out', 'linewidth', 2, 'fontsize', 10);

hold on

S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 0.2);
end

hold on
% Plotting Class 1 stations with custom colormap for stations
for i = 1:size(class1Stations, 1)
    colorIdx = interp1(linspace(min(class1Depths), max(class1Depths), size(stationCmap, 1)), 1:size(stationCmap, 1), class1Depths(i), 'nearest', 'extrap');
    depthColor = stationCmap(colorIdx, :);

     % Check the cluster number and set the marker style
    if class1Stations.Cluster(i) == 3
        markerStyle = 's'; % square for cluster 3
    else
        markerStyle = 'o'; % circle for clusters 1 and 2
    end

    m_plot(class1Stations.Longitude(i), class1Stations.Latitude(i), markerStyle, 'Color', depthColor, 'MarkerFaceColor', depthColor, 'MarkerSize', 8.5, 'MarkerEdgeColor', 'k'); 
end

% Create invisible axes for the second colorbar
secondAx = axes('Position', mainAx.Position, 'Visible', 'off', 'Color', 'none');
colormap(secondAx, stationCmap); % Set colormap for station depths
caxis(secondAx, [min(class1Depths) max(class1Depths)]);

cb2 = colorbar(secondAx, 'southoutside');
ylabel(cb2, 'Depth (km)');
cb2.FontSize = 13;
cb2.Position = [0.1, 0.33, 0.6, 0.03];

%linkaxes([mainAx, secondAx], 'xy');
%print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/figureS5a','-vector','-dpdf','-r0');
%%
% Plot Class 2 Stations Colored by Thickness
figure(2); clf;

% Create main axes for plotting
mainAx = axes('Position', [0.1, 0.3, 0.6, 0.6]);

m_scatter(velocityData.longitude, velocityData.latitude, 250, velocityData.vsv, 'Marker', '.');
colormap(customColormap2);
caxis([min(velocityData.vsv), max(velocityData.vsv)]); % Adjust as needed
cb1.Position = [0.75, 0.45, 0.03, 0.5];
cb1 = colorbar(mainAx,'eastoutside');

cb1.FontSize = 13;
freezeColors(mainAx);  % Freeze the colors of the plot
cbfreeze(cb1); % Freeze the colors of the first colorbar
ylabel(cb1, 'Velocity');

m_proj('Mercator', 'lon', [-128 -63], 'lat', [25 49.5]);
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir', 'out', 'linewidth', 2, 'fontsize', 10);

hold on

S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 0.2);
end

hold on;

% Plotting Class 2 stations with custom colormap
for i = 1:size(class2Stations, 1)
    colorIdx = interp1(linspace(min(class2Thicknesses), max(class2Thicknesses), size(stationCmap, 1)), 1:size(stationCmap, 1), class2Thicknesses(i), 'nearest', 'extrap');
    thicknessColor = stationCmap(colorIdx, :);
    m_plot(class2Stations.Longitude(i), class2Stations.Latitude(i), class2Shapes(i), 'Color', thicknessColor, 'MarkerFaceColor', thicknessColor, 'MarkerSize', 8.5, 'MarkerEdgeColor', 'k'); 
end

% Create invisible axes for the second colorbar
secondAx = axes('Position', mainAx.Position, 'Visible', 'off', 'Color', 'none');
colormap(secondAx, stationCmap); % Set colormap for station depths
caxis([min(class2Thicknesses) max(class2Thicknesses)]);

cb2 = colorbar(secondAx, 'southoutside');
ylabel(cb2, 'Thickness (km)');
cb2.FontSize = 13;
cb2.Position = [0.1, 0.33, 0.6, 0.03];
ylabel(cb2, 'Thickness (km)');

%print(figure(2),'/Users/evets/Desktop/Earthscope/Figures/figureS5b','-vector','-dpdf','-r0');
%%
% Calculate and display mean thickness value for Class 2
meanThickness = mean(class2Thicknesses);
%%
% Assume cluster3Depths contains the depths for cluster 3 stations
cluster3Depths = class1Depths(class1Stations.Cluster == 3);

% Find the starting depth for cluster 3, which is the minimum depth of cluster 3 stations
cluster3StartDepth = min(cluster3Depths);

% Create histograms for Class 1 depths and Class 2 thicknesses
figure(3); clf; % Figure for Class 1 depths histogram
histogram(class1Depths, 'BinMethod', 'auto', 'FaceColor', 'b','BinWidth',10);
xlabel('Depth (km)');
ylabel('Frequency');
title('Histogram of Class 1 Station Depths');

% Add a vertical dashed line for cluster 3 start depth
hold on; % Keep the histogram plot
line([cluster3StartDepth cluster3StartDepth], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Add legend to identify the line
legend('Class 1 Depths', 'Start of Cluster 3 Depths');
%%
squareThicknesses = class2Thicknesses(class2Shapes == 's');
circleThicknesses = class2Thicknesses(class2Shapes == 'o');

% Create histograms for "square" and "circle" stations
figure(4); clf;

% Histogram for Class 2 "square" stations
subplot(1, 2, 1);
histogram(squareThicknesses, 'FaceColor', 'b', 'BinWidth', 7);
xlabel('Thickness (km)');
ylabel('Frequency');
title('Class 2 "Square" Station Thicknesses');

% Histogram for Class 2 "circle" stations
subplot(1, 2, 2);
histogram(circleThicknesses, 'FaceColor', 'b', 'BinWidth', 7);
xlabel('Thickness (km)');
ylabel('Frequency');
title('Class 2 "Circle" Station Thicknesses');
%print(figure(4),'/Users/evets/Desktop/Earthscope/Figures/figureS5d','-vector','-dpdf','-r0');