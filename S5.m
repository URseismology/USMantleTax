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
% Plot Class 1 Stations Colored by Depth
figure(1); clf;
m_proj('Mercator', 'lon', [-128 -63], 'lat', [25 49.5]);
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir', 'out', 'linewidth', 2, 'fontsize', 10);

S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 0.2);
end

hold on;

% Plotting Class 1 stations with custom colormap
for i = 1:size(class1Stations, 1)
    colorIdx = interp1(linspace(min(class1Depths), max(class1Depths), size(customColormap, 1)), 1:size(customColormap, 1), class1Depths(i), 'nearest', 'extrap');
    depthColor = customColormap(colorIdx, :);
    m_plot(class1Stations.Longitude(i), class1Stations.Latitude(i), 'o', 'Color', depthColor, 'MarkerFaceColor', depthColor, 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
end

% Adjust colorbar for depth with custom colormap
colormap(customColormap);
cb=colorbar('southoutside');
cb.FontSize = 13;
caxis([min(class1Depths) max(class1Depths)]);
ylabel(cb, 'Depth (km)');
%print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/figureS5a','-vector','-dpdf','-r0');
%%
% Plot Class 2 Stations Colored by Thickness
figure(2); clf;
m_proj('Mercator', 'lon', [-128 -63], 'lat', [25 49.5]);
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir', 'out', 'linewidth', 2, 'fontsize', 10);

S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 0.2);
end
hold on;

% Plotting Class 2 stations with custom colormap
for i = 1:size(class2Stations, 1)
    colorIdx = interp1(linspace(min(class2Thicknesses), max(class2Thicknesses), size(customColormap, 1)), 1:size(customColormap, 1), class2Thicknesses(i), 'nearest', 'extrap');
    thicknessColor = customColormap(colorIdx, :);
    m_plot(class2Stations.Longitude(i), class2Stations.Latitude(i), class2Shapes(i), 'Color', thicknessColor, 'MarkerFaceColor', thicknessColor, 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
end

% Adjust colorbar for thickness with custom colormap
colormap(customColormap);
cb=colorbar('southoutside');
cb.FontSize = 13;
caxis([min(class2Thicknesses) max(class2Thicknesses)]);
ylabel(cb, 'Thickness (km)');

% Calculate and display mean thickness value for Class 2
meanThickness = mean(class2Thicknesses);
text('Position', [72, 35], 'String', ['Mean Thickness: ', num2str(meanThickness, '%.2f'), ' km'], 'FontSize', 12, 'Color', 'white', 'BackgroundColor', 'black', 'HorizontalAlignment', 'left');
%print(figure(2),'/Users/evets/Desktop/Earthscope/Figures/figureS5b','-vector','-dpdf','-r0');