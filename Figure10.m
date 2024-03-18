addpath(genpath(pwd));
clear; clc;

% Load the data
NVG_data = readtable('NVG_depth_data2.csv');
positives_data = readtable('positives_data.txt');

% Read provinces boundary coordinates from the CSV file
Basin_and_range = readtable('Basin_and_range.csv');
Appalatians = readtable('Appalatians.csv');
Snake_river_plain = readtable('Snake_river_plain.csv');
Colorado_plateau = readtable('Colorado_plateau.csv');
Rocky_mountains = readtable('Rocky_mountains.csv');

% Read US geological provinces
S = shaperead('physio.shp');

%Archean-proterozoic boundary
Archean_boundary = readtable('Archean_boundary.csv');

% Initialize arrays for storing class 1 and class 2 station data
SingleDiscontinuityStations = [];
SingleDiscontinuityDepths = [];
PairedLayerStations = [];
PairedLayerThicknesses = [];
PairedLayerShapes = []; % 'o' for circle, 's' for square
%%
% Classify and plot stations with single layers
for i = 1:height(NVG_data)
    station = NVG_data(i, :);
    station_name = station.StationName{1};
    station_cluster = station.Cluster;
    station_depth = station.Depth;

    idx_in_positives = find(strcmp(positives_data.StationName, station_name));

    % Retrieve cluster information for stations found in positives_data
    if ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
    else
        cluster_in_positives = [];
    end

    % Start conditions for single intra-lithosphere discontinuities with no base
    if isempty(idx_in_positives) && (station_cluster == 1 || station_cluster == 2 || station_cluster == 3)
        SingleDiscontinuityStations = [SingleDiscontinuityStations; station];
        SingleDiscontinuityDepths = [SingleDiscontinuityDepths; station_depth];
    elseif ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
        if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 3)
            SingleDiscontinuityStations = [SingleDiscontinuityStations; station];
            SingleDiscontinuityDepths = [SingleDiscontinuityDepths; station_depth];
        end
    end    
    %End conditions for single intra-lithosphere discontinuities with no base

    % % Start Condition for single transitional dsicontinuity with no base
    if station_cluster == 4 &&  ~any(ismember(cluster_in_positives, [1, 2]))
        SingleDiscontinuityStations = [SingleDiscontinuityStations; station];
        SingleDiscontinuityDepths = [SingleDiscontinuityDepths; station_depth];
    end
    %End Condition for single transitional layer
end
    

%%
% Classify and plot stations with layered structures
for i = 1:height(NVG_data)
    station = NVG_data(i, :);
    station_name = station.StationName{1};
    station_cluster = station.Cluster;
    station_depth = station.Depth;

    idx_in_positives = find(strcmp(positives_data.StationName, station_name));

    % Retrieve cluster information for stations found in positives_data
    if ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
    else
        cluster_in_positives = [];
    end

    % Start conditions for paired intra-lithospheric layers with top and bottom boundary
    if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 1 | cluster_in_positives == 2)
        for j = 1:length(idx_in_positives)
            depth_in_positives = positives_data.Depth(idx_in_positives);
            thickness = abs(station_depth - depth_in_positives);
            PairedLayerStations = [PairedLayerStations; station];
            PairedLayerThicknesses = [PairedLayerThicknesses; thickness];
            if station_depth < depth_in_positives(j)
                PairedLayerShapes = [PairedLayerShapes; 'o']; % P-type
            else
                PairedLayerShapes = [PairedLayerShapes; 's']; % N-type
            end
        end
    end
    %End conditions for paired intra-lithospheric discontinuities with top and bottom boundary

    %Start conditions for paired transitional discontinuities with top boundary
    if (station_cluster == 4) && any(cluster_in_positives == 1 | cluster_in_positives == 2) && ~any(cluster_in_positives == 3)
        PairedLayerStations = [PairedLayerStations; station];
        PairedLayerThicknesses = [PairedLayerThicknesses; thickness];
        PairedLayerShapes = [PairedLayerShapes; '^']; % Transitional with a top boundary
    end
    %End conditions for paired transitional discontinuities with top and bottom boundary

end
%%
%Creating custom colormap for depths, thicknesses and velocity
controlPoints = [0, 0, 0; ...   % black
                 0.5, 0, 0.5; ... % deep purple
                 1, 0.5, 0; ...   % orange
                 1, 1, 0.5];      % light yellow

% Create a linearly spaced vector for interpolation
x = linspace(1, size(controlPoints, 1), 256);  % 256 is the desired number of colors in the colormap

% Interpolate RGB values to create the colormap
DepthTicknessColormap = interp1(1:size(controlPoints, 1), controlPoints, x', 'linear');

% Custom colormap from mud brown to deep blue for velocities
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
VelocityColormap = interp1(1:size(controlPoints, 1), controlPoints, x, 'linear');

%stationCmap = customColormap;
stationCmap = flipud(DepthTicknessColormap);
%%
% Figure 10b
velocityData = readtable('Schultz_Pelkum_vpvs.csv');

% Set up figure and background
figure(1); clf;

% Create main axes for plotting
mainAx = axes('Position', [0.1, 0.3, 0.6, 0.6]);

% Plot seismic velocities
m_scatter(velocityData.longitude, velocityData.latitude, 250, velocityData.vsv, 'Marker', '.');
caxis([min(velocityData.vsv), max(velocityData.vsv)]); % Adjust as needed
colormap(mainAx,VelocityColormap);

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

% List of provinces to plot
provincesToPlot = {'COASTAL PLAIN','PACIFIC BORDER'};

% Loop through the shapes and plot only the specified provinces
for k = 1:length(S)
    if any(strcmp(S(k).PROVINCE, provincesToPlot))
        m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1.5);
    end
end

% % states
M=m_shaperead('ne_50m_admin_1_states_provinces');
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color', [0.6 0.6 0.6],'linewi',0.05);
end

m_line(Basin_and_range.Longitude, Basin_and_range.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Appalatians.Longitude, Appalatians.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Snake_river_plain.Longitude, Snake_river_plain.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Colorado_plateau.Longitude, Colorado_plateau.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Rocky_mountains.Longitude, Rocky_mountains.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% Overlay the archean proterozoic boundary as a red dashed line
m_line(Archean_boundary.Longitude, Archean_boundary.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 6);

hold on
%%
% Plotting stations with single discontinuities colored by depth
for i = 1:size(SingleDiscontinuityStations, 1)
    colorIdx = interp1(linspace(min(SingleDiscontinuityDepths), max(SingleDiscontinuityDepths), size(stationCmap, 1)), ...
        1:size(stationCmap, 1), SingleDiscontinuityDepths(i), 'nearest', 'extrap');
    depthColor = stationCmap(colorIdx, :);

     % Check the cluster number and set the marker style
    if SingleDiscontinuityStations.Cluster(i) == 4
        markerStyle = 's'; % square for cluster 3
    else
        markerStyle = 'o'; % circle for clusters 1 and 2
    end

    m_plot(SingleDiscontinuityStations.Longitude(i), SingleDiscontinuityStations.Latitude(i), markerStyle, 'Color', depthColor, ...
        'MarkerFaceColor', depthColor, 'MarkerSize', 8.5, 'MarkerEdgeColor', 'k'); 
end

% Create invisible axes for the second colorbar
secondAx = axes('Position', mainAx.Position, 'Visible', 'off', 'Color', 'none');
colormap(secondAx, stationCmap); % Set colormap for station depths
caxis(secondAx, [min(SingleDiscontinuityDepths) max(SingleDiscontinuityDepths)]);

cb2 = colorbar(secondAx, 'southoutside');
ylabel(cb2, 'Depth (km)');
cb2.FontSize = 13;
cb2.Position = [0.1, 0.33, 0.6, 0.03];

%linkaxes([mainAx, secondAx], 'xy');
print(figure(2),'/Users/evets/Desktop/Earthscope/Figures/figure10a','-vector','-dpdf','-r0');
%%
% Figure 10a
figure(2); clf;

% Create main axes for plotting
mainAx = axes('Position', [0.1, 0.3, 0.6, 0.6]);

m_scatter(velocityData.longitude, velocityData.latitude, 250, velocityData.vsv, 'Marker', '.');
colormap(VelocityColormap);
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

% Loop through the shapes and plot only the specified provinces
for k = 1:length(S)
    if any(strcmp(S(k).PROVINCE, provincesToPlot))
        m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1.5);
    end
end

% % states
M=m_shaperead('ne_50m_admin_1_states_provinces');
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color', [0.6 0.6 0.6],'linewi',0.05);
end

m_line(Basin_and_range.Longitude, Basin_and_range.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Appalatians.Longitude, Appalatians.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Snake_river_plain.Longitude, Snake_river_plain.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Colorado_plateau.Longitude, Colorado_plateau.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
m_line(Rocky_mountains.Longitude, Rocky_mountains.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% Overlay the archean proterozoic boundary as a red dashed line
m_line(Archean_boundary.Longitude, Archean_boundary.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 6);

hold on;

% Plot Stations with paired layers Colored by Thickness
for i = 1:size(PairedLayerStations, 1)
    colorIdx = interp1(linspace(min(PairedLayerThicknesses), max(PairedLayerThicknesses), ...
        size(stationCmap, 1)), 1:size(stationCmap, 1), PairedLayerThicknesses(i), 'nearest', 'extrap');
    thicknessColor = stationCmap(colorIdx, :);
    m_plot(PairedLayerStations.Longitude(i), PairedLayerStations.Latitude(i), ...
        PairedLayerShapes(i), 'Color', thicknessColor, 'MarkerFaceColor', thicknessColor, 'MarkerSize', 8.5, 'MarkerEdgeColor', 'k'); 
end

% Create invisible axes for the second colorbar
secondAx = axes('Position', mainAx.Position, 'Visible', 'off', 'Color', 'none');
colormap(secondAx, stationCmap); % Set colormap for station depths
caxis([min(PairedLayerThicknesses) max(PairedLayerThicknesses)]);

cb2 = colorbar(secondAx, 'southoutside');
ylabel(cb2, 'Thickness (km)');
cb2.FontSize = 13;
cb2.Position = [0.1, 0.33, 0.6, 0.03];
ylabel(cb2, 'Thickness (km)');

print(figure(4),'/Users/evets/Desktop/Earthscope/Figures/figure10a_inset','-vector','-dpdf','-r0');

%%
% Calculate and display mean thickness value for paired layers
meanThickness = mean(PairedLayerThicknesses);

% Assume cluster3Depths contains the depths for cluster 3 stations
cluster3Depths = SingleDiscontinuityDepths(SingleDiscontinuityStations.Cluster == 4);

% Find the starting depth for cluster 3, which is the minimum depth of cluster 3 stations
cluster4StartDepth = min(cluster3Depths);

%Figure10b inset
% Create histograms for Class 1 depths and Class 2 thicknesses
figure(3); clf; % Figure for Class 1 depths histogram
histogram(SingleDiscontinuityDepths, 'BinMethod', 'auto', 'FaceColor', 'b','BinWidth',10);
xlabel('Depth (km)');
ylabel('Frequency');
title('Histogram of Class 1 Station Depths');

% Add a vertical dashed line for cluster 3 start depth
hold on; % Keep the histogram plot
line([cluster4StartDepth cluster4StartDepth], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);

% Add legend to identify the line
legend('Class 1 Depths', 'Start of Cluster 4 Depths');
%%
squareThicknesses = PairedLayerThicknesses(PairedLayerShapes == 's');
circleThicknesses = PairedLayerThicknesses(PairedLayerShapes == 'o');

%%Figure10a inset
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