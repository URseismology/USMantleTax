clear;clc;
addpath(genpath(pwd));

% Read provinces boundary coordinates from the CSV file
Basin_and_range = readtable('Basin_and_range.csv');
Appalatians = readtable('Appalatians.csv');
Snake_river_plain = readtable('Snake_river_plain.csv');
Colorado_plateau = readtable('Colorado_plateau.csv');
Rocky_mountains = readtable('Rocky_mountains.csv');

%Archean-proterozoic boundary
Archean_boundary = readtable('Archean_boundary.csv');
%%
% Initialize base map
figure(1); clf;

hh = figure (1); clf;
fs=10;
m_proj('Mercator','lon',[-126 -63],'lat',[25 49.5]); %continental US
shading flat;

% Read US geological provinces
S = shaperead('physio.shp');

% List of provinces to plot
provincesToPlot = {'COASTAL PLAIN','PACIFIC BORDER'};

% Loop through the shapes and plot only the specified provinces
for k = 1:length(S)
    if any(strcmp(S(k).PROVINCE, provincesToPlot))
        m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 3);
    end
end

m_line(Basin_and_range.Longitude, Basin_and_range.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 3);
m_line(Appalatians.Longitude, Appalatians.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 3);
m_line(Snake_river_plain.Longitude, Snake_river_plain.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 3);
m_line(Colorado_plateau.Longitude, Colorado_plateau.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 3);
m_line(Rocky_mountains.Longitude, Rocky_mountains.Latitude, 'color', 'k', 'LineStyle', '-', 'LineWidth', 3);

% Overlay the archean proterozoic boundary as a red dashed line
m_line(Archean_boundary.Longitude, Archean_boundary.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 6);

% Coastlines and Grids
m_coast('line','color','k','linewidth',2);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% % states
M=m_shaperead('ne_50m_admin_1_states_provinces');
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color', [0.6 0.6 0.6],'linewi',0.05);
end
%%
hold on
% Load the data
NVG_data = readtable('NVG_depth_data2.csv');
positives_data = readtable('positives_data.txt');

% Initialize counters for each class and case
Intra_with_no_base_inside = 0; Intra_with_no_base_outside = 0;
Intra_with_no_base_total = 0;
Transitional_with_no_top_inside = 0; Transitional_with_no_top_outside = 0;
Transitional_with_no_top_total = 0;
P_type_layer_inside = 0; P_type_layer_outside = 0;
P_type_layer_total = 0;
N_type_layer_inside = 0; N_type_layer_outside = 0;
N_type_layer_total = 0;
transitional_with_top_inside = 0; transitional_with_top_outside = 0;
transitional_with_top_total = 0;
sublithospheric_inside = 0; sublithospheric_outside = 0;
sublithospheric_total = 0;

Intra_with_no_base = 0; %No top or bottom boundary
N_type_layer = 0; % NVG depth < positives depth
P_type_layer = 0; % NVG depth > positives depth
X_discontinuity = 0; % Counter for Cluster 3 and 4 transitional stations with no boundary
Transitional_discontinuities_with_top = 0; % Counter for Cluster 4 discontinuities with a specific condition


% Classify and plot stations
for i = 1:height(NVG_data)
    station = NVG_data(i, :);
    station_name = station.StationName{1};
    station_cluster = station.Cluster;
    station_depth = station.Depth;

    % Determine if the station is inside the Archean-Proterozoic boundary
    insideBoundary = isInsideBoundary(station.Longitude, station.Latitude, Archean_boundary.Longitude, Archean_boundary.Latitude);

    idx_in_positives = find(strcmp(positives_data.StationName, station_name));

    % Retrieve cluster information for stations found in positives_data
    if ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
    else
        cluster_in_positives = [];
    end

    % Start conditions for single intra-lithosphere discontinuities with no base
    if isempty(idx_in_positives) && (station_cluster == 1 || station_cluster == 2 || station_cluster == 3)
        if insideBoundary
            m_plot(station.Longitude, station.Latitude, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            Intra_with_no_base_inside = Intra_with_no_base_inside + 1;
        else
            m_plot(station.Longitude, station.Latitude, 'o', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            Intra_with_no_base_outside = Intra_with_no_base_outside + 1;
        end
        Intra_with_no_base_total = Intra_with_no_base_total + 1;

    elseif ~isempty(idx_in_positives)
        cluster_in_positives = positives_data.Cluster(idx_in_positives);
        if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 3)
            if insideBoundary
                m_plot(station.Longitude, station.Latitude, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                Intra_with_no_base_inside = Intra_with_no_base_inside + 1;
            else
                m_plot(station.Longitude, station.Latitude, 'o', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                Intra_with_no_base_outside = Intra_with_no_base_outside + 1;
            end

            Intra_with_no_base_total = Intra_with_no_base_total + 1;
        end
    end
    %End conditions for single intra-lithosphere discontinuities with no base

    % Start Condition for single transitional dsicontinuity with no base
    if station_cluster == 4 &&  ~any(ismember(cluster_in_positives, [1, 2]))

        % Choose plotting color and marker based on whether it's inside or outside the boundary
        if insideBoundary
            plotColor = 'y'; % Replace 'uniqueColor1' with the desired color for inside boundary
            plotMarker = 'o'; % Example marker
            Transitional_with_no_top_inside = Transitional_with_no_top_inside + 1;
        else
            plotColor = [0.65 0.4 0.2]; % Replace 'uniqueColor2' with the desired color for outside boundary
            plotMarker = 'o'; % Example marker
            Transitional_with_no_top_outside = Transitional_with_no_top_outside + 1;
        end
        m_plot(station.Longitude, station.Latitude, plotMarker, 'MarkerFaceColor', plotColor, 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
        % Increment a specific counter if needed

        Transitional_with_no_top_total = Transitional_with_no_top_total+1;
    end
    %End Condition for single transitional layer

    % Start conditions for paired intra-lithospheric layers with top and bottom boundary
    if (station_cluster == 1 || station_cluster == 2 || station_cluster == 3) && any(cluster_in_positives == 1 | cluster_in_positives == 2)
        for j = 1:length(idx_in_positives)
            depth_in_positives = positives_data.Depth(idx_in_positives);
            if station_depth < depth_in_positives(j)
                if insideBoundary
                    m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', 'g', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                    P_type_layer_inside = P_type_layer_inside + 1;
                else
                    m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                    P_type_layer_outside = P_type_layer_outside + 1;
                end
                P_type_layer_total = P_type_layer_total + 1;
                break; % Prevent plotting the same station multiple times
            else
                if insideBoundary
                    m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', 'g', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                    N_type_layer_inside = N_type_layer_inside + 1;
                else
                    m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
                    N_type_layer_outside = N_type_layer_outside + 1;
                end
                N_type_layer_total = N_type_layer_total + 1;
                break; % Prevent plotting the same station multiple times
            end
        end
        %End conditions for paired intra-lithospheric discontinuities with top and bottom boundary
    end

    %Start conditions for paired transitional discontinuities with top boundary
    if (station_cluster == 4) && any(cluster_in_positives == 1 | cluster_in_positives == 2) && ~any(cluster_in_positives == 3)
        if insideBoundary
            m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', 'y', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            transitional_with_top_inside = transitional_with_top_inside + 1;
        else
            m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            transitional_with_top_outside = transitional_with_top_outside + 1;
        end

        transitional_with_top_total = transitional_with_top_total + 1;
    end
    %End conditions for paired transitional discontinuities with top and bottom boundary
end

% Single Sub_lithospheric discontinuities
for i = 1:height(positives_data)
    station = positives_data(i, :);
    insideBoundary = isInsideBoundary(station.Longitude, station.Latitude, Archean_boundary.Longitude, Archean_boundary.Latitude);
    if positives_data.Cluster(i) == 4
        positive_station_name = positives_data.StationName{i};
        if insideBoundary
            m_plot(positives_data.Longitude(i), positives_data.Latitude(i), 'o', 'MarkerFaceColor', 'magenta', 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            sublithospheric_inside = sublithospheric_inside + 1;
        else
            m_plot(station.Longitude, station.Latitude, 's', 'MarkerFaceColor', [0.65 0.4 0.2], 'MarkerSize', 11, 'MarkerEdgeColor', 'k');
            sublithospheric_outside = sublithospheric_outside + 1;
        end
        sublithospheric_total = sublithospheric_total+ 1;
    end
end

%%
% Output the counts
fprintf('Intra_with_no_base_inside: %d\n', Intra_with_no_base_inside);
fprintf('Intra_with_no_base_outside: %d\n', Intra_with_no_base_outside );
fprintf('Intra_with_no_base_total: %d\n', Intra_with_no_base_total);
fprintf('Transitional_with_no_top_inside: %d\n', Transitional_with_no_top_inside);
fprintf('Transitional_with_no_top_outside: %d\n', Transitional_with_no_top_outside);
fprintf('Transitional_with_no_top_total: %d\n', Transitional_with_no_top_total);
fprintf('P_type_layer_inside: %d\n', P_type_layer_inside);
fprintf('P_type_layer_outside: %d\n', P_type_layer_outside);
fprintf('P_type_layer_total: %d\n', P_type_layer_total);
fprintf('N_type_layer_inside: %d\n', N_type_layer_inside);
fprintf('N_type_layer_outside: %d\n', N_type_layer_outside);
fprintf('N_type_layer_total: %d\n', N_type_layer_total);
fprintf('transitional_with_top_inside: %d\n', transitional_with_top_inside);
fprintf('transitional_with_top_outside: %d\n', transitional_with_top_outside);
fprintf('transitional_with_top_total: %d\n', transitional_with_top_total);
fprintf('sublithospheric_inside: %d\n', sublithospheric_inside);
fprintf('sublithospheric_outside: %d\n', sublithospheric_outside);
fprintf('sublithospheric_total: %d\n', sublithospheric_total);


hold off

%print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/figure9','-vector','-dpdf','-r0');

