addpath(genpath(pwd));
clear;clc;

% Load the data
NVG_data = readtable('NVG_depth_data2.csv');
positives_data = readtable('positives_data.txt');

% Initialize counters for each class and case
count_class_1 = 0;
count_class_2_lower = 0; % NVG depth < positives depth
count_class_2_higher = 0; % NVG depth > positives depth
count_class_3_NVG = 0;
count_class_3_positives = 0;

% Initialize map
figure(1); clf;
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]); %continental US
shading flat;

% Read and plot the physiographic provinces shapefile
S = shaperead('physio.shp');
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Coastlines and Grids
m_coast('line', 'color', 'k', 'linewidth', 1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

hold on
% Classify and plot stations
for i = 1:height(NVG_data)
    station = NVG_data(i, :);
    station_name = station.StationName{1};
    station_cluster = station.Cluster;
    station_depth = station.Depth;

    idx_in_positives = find(strcmp(positives_data.StationName, station_name));

    % Class 3 stations from NVG data
    if station_cluster == 4
        m_plot(station.Longitude, station.Latitude, 'o', 'Color', [1 0.5 0], 'MarkerFaceColor', [1 0.5 0], 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
        count_class_3_NVG = count_class_3_NVG + 1;
    elseif isempty(idx_in_positives) && (station_cluster == 1 || station_cluster == 2)
        % Class 1
        m_plot(station.Longitude, station.Latitude, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
        count_class_1 = count_class_1 + 1;
    elseif ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
        depth_in_positives = positives_data.Depth(idx_in_positives);

        % Class 1 and Class 2 logic
        if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 3)
            % Class 1
            m_plot(station.Longitude, station.Latitude, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
            count_class_1 = count_class_1 + 1;
        elseif (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 1 | cluster_in_positives == 2)
            % Class 2
            for j = 1:length(idx_in_positives)
                if station_depth < depth_in_positives(j)
                    m_plot(station.Longitude, station.Latitude, 'o', 'MarkerFaceColor', 'y', 'MarkerSize', 9, 'MarkerEdgeColor', 'k'); 
                    count_class_2_lower = count_class_2_lower + 1;
                    break; % Prevent plotting the same station multiple times
                else
                    m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', 'y', 'MarkerSize', 9, 'MarkerEdgeColor', 'k');
                    count_class_2_higher = count_class_2_higher + 1;
                    break; % Prevent plotting the same station multiple times
                end
            end
        end
    end
end

% Classify and plot Class 3 stations from positives_data
for i = 1:height(positives_data)
    if positives_data.Cluster(i) == 4
        m_plot(positives_data.Longitude(i), positives_data.Latitude(i), 'o', 'Color', [1 0 1], 'MarkerFaceColor', [1 0 1], 'MarkerSize', 9, 'MarkerEdgeColor', 'k');
        count_class_3_positives = count_class_3_positives + 1;
    end
end

% Output the counts
fprintf('Number of Class 1 stations: %d\n', count_class_1);
fprintf('Number of Class 2 stations (Lower Depth): %d\n', count_class_2_lower);
fprintf('Number of Class 2 stations (Higher Depth): %d\n', count_class_2_higher);
fprintf('Number of Class 3 stations (NVG data): %d\n', count_class_3_NVG);
fprintf('Number of Class 3 stations (Positives data): %d\n', count_class_3_positives);

hold off

%print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/figure10','-vector','-dpdf','-r0');

